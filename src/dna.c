#include "dna.h"

PG_MODULE_MAGIC;


/**
 * @brief Structure used to store the state of the K-mer generator.
 */
typedef struct KmerGeneratorState {
    uint8_t kmer_length;       /**< Length of the K-mer(s) to generate */
    uint32_t length;           /**< Total length of the DNA sequence */
    uint8_t* byte_ptr;         /**< Pointer to a byte in the byte array representing the DNA sequence */
    uint8_t nucleotide_ctr;    /**< Counter for nucleotides in the current byte (0-3) */
} KmerGeneratorState;


/**
 * @brief Adds the length of the last byte of the DNA sequence to the beginning of the DNA object.
 * This is necessary because the last byte might not hold 4 nucleotides. 
 * 
 * @param dna The DNA object.
 * @param length The amount of nucleotides of the DNA sequence.
 */
static void store_last_byte_length(DNA* dna, uint32_t length) {
    uint8_t* data_ptr = (uint8_t*) VARDATA(dna); 
    *data_ptr = length % 4 == 0 ? 0b00000100 : length % 4;
}

/**
 * @brief Creates a DNA object from a string.
 * 
 * @param str The string representing the DNA sequence.
 * @param length The length of the DNA sequence.
 * @return A pointer to the created DNA object.
 */
static DNA* make_dna(const char* str, uint32_t length) {
    DNA* dna = palloc0(VARHDRSZ + ceil((double) length / 4) + 1);     // + 1 for last byte length
    SET_VARSIZE(dna, VARHDRSZ + ceil((double) length / 4) + 1);
    uint8_t* data_ptr = (uint8_t*) VARDATA(dna);

    store_last_byte_length(dna, length);
    data_ptr++;                                                       // skip first byte for last byte length

    for (uint32_t i = 0; i < length; i += 4) {
        uint8_t current_byte = 0b00000000;                            // initialize byte with 0

        for (char j = 0; j < 4; j++) {                                // read 4 nucleotides at a time
            char c = str[i + j];
            if (i + j >= length) {                                    // shift byte to the left by 2 bits for each missing nucleotide
                current_byte <<= 2;                                   // This is necessary because the last byte might not hold 4 nucleotides  
                continue;
            }
            switch (toupper(c)) {
                case 'A': current_byte = (current_byte << 2) | 0b00; break; // TODO: Use LUT for this
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

/**
 * @brief Parses a DNA sequence from a string.
 * 
 * @param str The string representing the DNA sequence to parse.
 * @return A pointer to the DNA object created from the string
 */
static DNA* dna_parse(const char* str) {
    uint32_t length = strlen(str);
    if (length == 0) {
		ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
	  	errmsg("kmer should not be empty")));
	}
    return make_dna(str, length);
}

/**
 * @brief Converts a DNA object to a string representation.
 * 
 * @param dna The DNA object to convert.
 * @return A string representation of the DNA sequence.
 */
static char* dna_to_string(DNA* dna) {
    // length represents the max amount of nucleotides in the DNA sequence, but the last byte might not be full
    uint32_t length = (VARSIZE(dna) - VARHDRSZ - 1) * 4;                               // - 1 for last byte length

    uint8_t* data_ptr = (uint8_t*) VARDATA(dna);
    uint8_t last_byte_length = *data_ptr++;
    char* str = palloc0(length + 1 - (4 - last_byte_length));

    for (uint32_t i = 0; i < length; i += 4) {
        uint8_t current_byte = *data_ptr;
        data_ptr++;
        uint8_t number_of_nucleotides_to_read = i + 4 == length ? last_byte_length : 4; // read less nucleotides for the last byte

        for (char j = 0; j < number_of_nucleotides_to_read; j++) {
            uint8_t nucleotide = (current_byte >> (6 - j * 2)) & 0b11;                  // shift to the right by 6 bits for the first nucleotide, 4 bits for the second, etc.
            if (i + j >= length) {                                                      // stop reading nucleotides if we reached the end of the DNA sequence
                break;
            }
            switch (nucleotide) {
                case 0b00: str[i + j] = 'A'; break; // TODO: Use LUT for this
                case 0b01: str[i + j] = 'C'; break;
                case 0b10: str[i + j] = 'G'; break;
                case 0b11: str[i + j] = 'T'; break;
            }
        }
    }
    str[length] = '\0';
    return str;
}


/**
 * @brief Gets the length of the DNA sequence.
 * 
 * @param dna The DNA object.
 * @return The length of the DNA sequence.
 */
static uint32_t get_dna_sequence_length(DNA* dna) {
    uint32_t length_in_bytes = (VARSIZE(dna) - VARHDRSZ - 1) * 4;
    uint8_t* data_ptr = (uint8_t*) VARDATA(dna);
    uint8_t last_byte_length = *data_ptr;
    return length_in_bytes - (4 - last_byte_length);
}

/**
 * @brief Generates a K-mer from a subsequence of the DNA sequence.
 * 
 * @param byte_ptr Pointer to the byte in the DNA sequence where the subsequence starts.
 * @param nucleotide_ctr Counter for nucleotides in the current byte (0-3) where the subsequence starts.
 * @param kmer_length The length of the K-mer to generate.
 * @return A pointer to the generated K-mer. 
 */
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
        uint8_t nucleotide = current_byte >> (6 - nucleotide_ctr * 2) & 0b11;   // shift to the right by 6 bits for the first nucleotide, 4 bits for the second, etc.
        kmer->value = kmer->value << 2 | nucleotide;                            // shift to the left by 2 bits for each nucleotide
        nucleotide_ctr++;
    }
    return kmer;
}

/**
 * @brief Initializes the state for K-mer generation.
 * 
 * @param state The KmerGeneratorState object to initialize.
 * @param dna The DNA object.
 * @param funcctx The function call context.
 */
static void init_kmer_generator_state(KmerGeneratorState *state, DNA *dna, FuncCallContext *funcctx) {
    if (state->kmer_length > 32) {
        errmsg("You cannot generate a kmer that is longer than 32 nucleotides.");
    }
    state->length = get_dna_sequence_length(dna);
    state->byte_ptr = (uint8_t *)VARDATA(dna);
    state->nucleotide_ctr = 0;
    state->byte_ptr++; // skip first byte for last byte length

    funcctx->max_calls = state->length - state->kmer_length + 1; // number of kmers to generate
}
/* ------------------------------------------------------------------------- */

/**
 * @brief Postgres input function for DNA.
 * 
 * @param str The input string.
 * @return The DNA object created from the input string.
 */
PG_FUNCTION_INFO_V1(dna_in);
Datum dna_in(PG_FUNCTION_ARGS) {
    char* str = PG_GETARG_CSTRING(0);
    DNA* result = dna_parse(str);
    PG_FREE_IF_COPY(str, 0);
    PG_RETURN_BYTEA_P(result);
}

/**
 * @brief Postgres output function for DNA.
 * 
 * @param dna The DNA object.
 * @return The string representation of the DNA sequence.
 */
PG_FUNCTION_INFO_V1(dna_out);
Datum dna_out(PG_FUNCTION_ARGS) {
    DNA* dna = PG_GETARG_BYTEA_P(0);
    //! elog(INFO, "dna_out: %x", dna);
    char* str = dna_to_string(dna);
    PG_FREE_IF_COPY(dna, 0);
    PG_RETURN_CSTRING(str);
}

/**
 * @brief Postgres send function for DNA.
 * 
 * @param dna The DNA object.
 * @return The bytea representation of the DNA object.
 */
PG_FUNCTION_INFO_V1(dna_recv);
Datum dna_recv(PG_FUNCTION_ARGS) {
    StringInfo buf = (StringInfo) PG_GETARG_POINTER(0);
    int32 len = pq_getmsgint(buf, 4);
    bytea* result = (bytea*) palloc(len + VARHDRSZ);
    SET_VARSIZE(result, len + VARHDRSZ);
    memcpy(VARDATA(result), pq_getmsgbytes(buf, len), len);
    PG_RETURN_BYTEA_P(result);
}

/**
 * @brief Postgres receive function for DNA.
 * 
 * @param dna The DNA object.
 * @return The bytea representation of the DNA object.
 */
PG_FUNCTION_INFO_V1(dna_send);
Datum dna_send(PG_FUNCTION_ARGS) {
    bytea* dna = PG_GETARG_BYTEA_P(0);
    StringInfoData buf;
    pq_begintypsend(&buf);
    pq_sendbytes(&buf, VARDATA(dna), VARSIZE_ANY_EXHDR(dna));
    PG_FREE_IF_COPY(dna, 0);
    PG_RETURN_BYTEA_P(pq_endtypsend(&buf));
}

/**
 * @brief Postgres cast function from text to DNA.
 * 
 * @param txt The text to cast.
 * @return The DNA object created from the text.
 */
PG_FUNCTION_INFO_V1(DNA_cast_from_text);
Datum DNA_cast_from_text(PG_FUNCTION_ARGS) {
    text *txt = PG_GETARG_TEXT_P(0);
    char *str = DatumGetCString(DirectFunctionCall1(textout,
               PointerGetDatum(txt)));
    DNA* dna = dna_parse(str);
    PG_FREE_IF_COPY(txt, 0);
    PG_RETURN_BYTEA_P(dna);
}

/**
 * @brief Postgres cast function from DNA to text.
 * 
 * @param dna The DNA object to cast.
 * @return The text representation of the DNA sequence.
 */
PG_FUNCTION_INFO_V1(DNA_cast_to_text);
Datum DNA_cast_to_text(PG_FUNCTION_ARGS) {
    DNA* dna = PG_GETARG_BYTEA_P(0);
    text* out = (text *)DirectFunctionCall1(textin,
              PointerGetDatum(dna_to_string(dna)));
    PG_FREE_IF_COPY(dna, 0);
    PG_RETURN_TEXT_P(out);
}


/**
 * @brief Postgres function to get the length of a DNA sequence.
 * 
 * @param dna The DNA object.
 * @return The length of the DNA sequence.
 */
PG_FUNCTION_INFO_V1(dna_length);
Datum dna_length(PG_FUNCTION_ARGS) {
    DNA* dna = PG_GETARG_BYTEA_P(0);
    uint32_t length = get_dna_sequence_length(dna);
    PG_FREE_IF_COPY(dna, 0);
    PG_RETURN_UINT32(length);
}

/**
 * @brief Postgres function to generate K-mers from a DNA sequence.
 * 
 * @param dna The DNA object.
 * @param kmer_length The length of the K-mers to generate.
 * @return A set of K-mers.
 */
PG_FUNCTION_INFO_V1(dna_generate_kmers);
Datum dna_generate_kmers(PG_FUNCTION_ARGS) {
    FuncCallContext *funcctx;

    if (SRF_IS_FIRSTCALL()) {
        MemoryContext oldcontext;
        funcctx = SRF_FIRSTCALL_INIT();
        oldcontext = MemoryContextSwitchTo(funcctx->multi_call_memory_ctx);

        // Store state that needs to persist across calls
        funcctx->user_fctx = palloc0(sizeof(KmerGeneratorState));
        KmerGeneratorState *state = (KmerGeneratorState *) funcctx->user_fctx;

        DNA* dna = PG_GETARG_BYTEA_P(0);
        state->kmer_length = (uint8_t) PG_GETARG_UINT16(1);
        init_kmer_generator_state(state, dna, funcctx);

        MemoryContextSwitchTo(oldcontext);
    }

    funcctx = SRF_PERCALL_SETUP();
    KmerGeneratorState *state = (KmerGeneratorState *) funcctx->user_fctx;

   
    if (funcctx->call_cntr < funcctx->max_calls) {
        Kmer* kmer = generate_kmer_subsequence(state->byte_ptr, state->nucleotide_ctr, state->kmer_length);
        
        if (++state->nucleotide_ctr == 4) {
            state->nucleotide_ctr = 0;
            state->byte_ptr++;
        }
        
        // Convert kmer to text for return
        SRF_RETURN_NEXT(funcctx, KmerPGetDatum(kmer));
    } else {
        SRF_RETURN_DONE(funcctx);
    }
}


