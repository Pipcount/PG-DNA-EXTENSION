#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "postgres.h"
#include "fmgr.h"

#include "libpq/pqformat.h"
#include "utils/fmgrprotos.h"
#include "utils/elog.h"
#include "utils/varlena.h"
#include "varatt.h"

PG_MODULE_MAGIC;

typedef bytea DNA;

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

static bytea* dna_parse(const char* str) {
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


PG_FUNCTION_INFO_V1(dna_in);
Datum dna_in(PG_FUNCTION_ARGS) {
    char* str = PG_GETARG_CSTRING(0);
    PG_RETURN_BYTEA_P(dna_parse(str));
}

PG_FUNCTION_INFO_V1(dna_out);
Datum dna_out(PG_FUNCTION_ARGS) {
    DNA* dna = PG_GETARG_BYTEA_P(0);
    elog(INFO, "dna_out: %x", dna);
    char* str = dna_to_string(dna);
    PG_FREE_IF_COPY(dna, 0);
    PG_RETURN_CSTRING(str);
}