#include "dna.h"

PG_MODULE_MAGIC;


static void store_last_byte_length(DNA* dna, uint32_t length) {
    uint8_t* data_ptr = (uint8_t*) VARDATA(dna); 
    *data_ptr = length % 4 == 0 ? 0b00000100 : length % 4;
}

static DNA* make_dna(const char* str, uint32_t length) {
    DNA* dna = palloc0(VARHDRSZ + ceil((double) length / 4) + 1);       // + 1 for last byte length
    SET_VARSIZE(dna, VARHDRSZ + ceil((double) length / 4) + 1);
    uint8_t* data_ptr = (uint8_t*) VARDATA(dna);
    store_last_byte_length(dna, length);
    data_ptr++;                                                                 // skip first byte for last byte length
    for (uint32_t i = 0; i < length; i += 4) {
        uint8_t current_byte = 0b00000000;

        for (char j = 0; j < 4; j++) {
            char c = str[i + j];
            if (i + j >= length) {
                current_byte <<= 2;
                continue;
            }
            switch (c) {
                case 'A': current_byte = (current_byte << 2) | 0b00; break;
                case 'C': current_byte = (current_byte << 2) | 0b01; break;
                case 'G': current_byte = (current_byte << 2) | 0b10; break;
                case 'T': current_byte = (current_byte << 2) | 0b11; break;
                default:
                    ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                    errmsg("invalid nucleotide")));
            }
        }
        *data_ptr = current_byte;
        data_ptr++;
    }
    return dna; 
}

static DNA* dna_parse(const char* str) {
    uint32_t length = strlen(str);
    if (length == 0) {
		ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
	  	errmsg("kmer should not be empty")));
	}
    return make_dna(str, length);
}

static char* dna_to_string(DNA* dna) {
    // length represents the max amount of nucleotides in the DNA sequence, but the last byte might not be full
    uint32_t length = (VARSIZE(dna) - VARHDRSZ - 1) * 4;        // - 1 for last byte length

    uint8_t* data_ptr = (uint8_t*) VARDATA(dna);
    uint8_t last_byte_length = *data_ptr++;
    char* str = palloc0(length + 1 - (4 - last_byte_length));

    for (uint32_t i = 0; i < length; i += 4) {
        uint8_t current_byte = *data_ptr;
        data_ptr++;
        uint8_t number_of_nucleotides_to_read = i + 4 == length ? last_byte_length : 4;
        for (char j = 0; j < number_of_nucleotides_to_read; j++) {
            uint8_t nucleotide = (current_byte >> (6 - j * 2)) & 0b11;
            if (i + j >= length) {
                break;
            }
            switch (nucleotide) {
                case 0b00: str[i + j] = 'A'; break;
                case 0b01: str[i + j] = 'C'; break;
                case 0b10: str[i + j] = 'G'; break;
                case 0b11: str[i + j] = 'T'; break;
            }
        }
    }
    str[length] = '\0';
    return str;
}

static uint32_t get_dna_sequence_length(DNA* dna) {
    uint32_t length_in_bytes = (VARSIZE(dna) - VARHDRSZ - 1) * 4;
    uint8_t* data_ptr = (uint8_t*) VARDATA(dna);
    uint8_t last_byte_length = *data_ptr;
    return length_in_bytes - (4 - last_byte_length);
}

static char* kmer_value_to_string(Kmer* kmer) {
	char str[32];
	//! elog(INFO, "kmer value (kmer_value_to_string): %lu", kmer -> value);
	for (uint8_t i = 0; i < kmer -> k; i++) {
        uint8_t shift = (kmer -> k - i - 1) * 2;
        uint8_t nucleotide = (kmer -> value >> shift) & 0b11;

        if (nucleotide == 0b00) {
            str[i] = 'A';
        } else if (nucleotide == 0b01) {
            str[i] = 'C';
        } else if (nucleotide == 0b10) {
            str[i] = 'G';
        } else if (nucleotide == 0b11) {
            str[i] = 'T';
        }
    }
	str[kmer->k] = '\0';
	return psprintf("%s", str);
}

static Kmer* generate_kmer_subsequence(uint8_t* byte_ptr, uint8_t nucleotide_ctr, uint8_t kmer_length) {
    Kmer* kmer = palloc0(sizeof(Kmer));
    kmer->value = 0;
    kmer->k = kmer_length;

    for (int i = 0; i < kmer_length; i++) {
        if (nucleotide_ctr == 4) {
            nucleotide_ctr = 0;
            (byte_ptr)++;
        }
        uint8_t current_byte = *byte_ptr;
        uint8_t nucleotide = current_byte >> (6 - nucleotide_ctr * 2) & 0b11;
        kmer->value = kmer->value << 2 | nucleotide;
        nucleotide_ctr++;
    }
    return kmer;
}


PG_FUNCTION_INFO_V1(dna_in);
Datum dna_in(PG_FUNCTION_ARGS) {
    char* str = PG_GETARG_CSTRING(0);
    PG_RETURN_BYTEA_P(dna_parse(str));
}

PG_FUNCTION_INFO_V1(dna_out);
Datum dna_out(PG_FUNCTION_ARGS) {
    DNA* dna = PG_GETARG_BYTEA_P(0);
    //! elog(INFO, "dna_out: %x", dna);
    char* str = dna_to_string(dna);
    PG_FREE_IF_COPY(dna, 0);
    PG_RETURN_CSTRING(str);
}

PG_FUNCTION_INFO_V1(dna_recv);
Datum dna_recv(PG_FUNCTION_ARGS) {
    StringInfo buf = (StringInfo) PG_GETARG_POINTER(0);
    int32 len = pq_getmsgint(buf, 4);
    bytea* result = (bytea*) palloc(len + VARHDRSZ);
    SET_VARSIZE(result, len + VARHDRSZ);
    memcpy(VARDATA(result), pq_getmsgbytes(buf, len), len);
    PG_RETURN_BYTEA_P(result);
}

PG_FUNCTION_INFO_V1(dna_send);
Datum dna_send(PG_FUNCTION_ARGS) {
    bytea* dna = PG_GETARG_BYTEA_P(0);
    StringInfoData buf;
    pq_begintypsend(&buf);
    pq_sendbytes(&buf, VARDATA(dna), VARSIZE_ANY_EXHDR(dna));
    PG_FREE_IF_COPY(dna, 0);
    PG_RETURN_BYTEA_P(pq_endtypsend(&buf));
}

PG_FUNCTION_INFO_V1(dna_length);
Datum dna_length(PG_FUNCTION_ARGS) {
    DNA* dna = PG_GETARG_BYTEA_P(0);
    uint32_t length = get_dna_sequence_length(dna);
    PG_FREE_IF_COPY(dna, 0);
    PG_RETURN_UINT32(length);
}

PG_FUNCTION_INFO_V1(dna_generate_kmers);
Datum dna_generate_kmers(PG_FUNCTION_ARGS) {
    // Declarations
    DNA* dna = PG_GETARG_BYTEA_P(0);
    uint8_t kmer_length = (uint8_t) PG_GETARG_UINT16(1);
    if (kmer_length > 32) {
        errmsg("You cannot generate a kmer that is longer than 32 nucleotides.");
    }
    uint32_t length = get_dna_sequence_length(dna);
    uint8_t* byte_ptr = (uint8_t*) VARDATA(dna);
    uint8_t nucleotide_ctr = 0;
    byte_ptr++; // skip first byte for last byte length

    for (uint32_t nb_kmers = 0; nb_kmers < (length - kmer_length + 1); nb_kmers++) {
        Kmer* kmer_subsequence = generate_kmer_subsequence(byte_ptr, nucleotide_ctr, kmer_length);
        if (++nucleotide_ctr == 4) {
            nucleotide_ctr = 0;
            byte_ptr++;
        }
        elog(INFO, "kmer: %s", kmer_value_to_string(kmer_subsequence));
    }

}
