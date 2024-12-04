#include "kmer.h"


/**
 * @brief Creates a K-mer from a string.
 * 
 * @param str The string representing the K-mer.
 * @param length The length of the K-mer.
 * @return A pointer to the created K-mer.
*/
static Kmer* make_kmer(const char *str, uint8_t length) {
	Kmer* kmer = palloc0(sizeof(Kmer));
	kmer -> k = length;
	kmer -> value = 0;

	for (uint8_t i = 0; i < length; i++) {
		char c = str[i];
		add_nucleotide_to_uint(kmer -> value, c);
	}
	return kmer;
}

/**
 * @brief Parses a K-mer from a string.
 * 
 * @param str The string representing the K-mer to parse.
 * @return A pointer to the K-mer created from the string.
 */
static Kmer* kmer_parse(const char* str) {
	uint8_t length = strlen(str);
	//! elog(INFO, "kmer length (kmer_parse): %d", length);
	if (length > 32) {
		ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
      	errmsg("kmer should not exceed 32 nucleotides")));
	} else if (length == 0) {
		ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
	  	errmsg("kmer should not be empty")));
	}
	return make_kmer(str, length);
}

/**
 * @brief Converts a K-mer to a string.
 * 
 * @param kmer The K-mer to convert.
 * @return The string representation of the K-mer.
 */
char* kmer_value_to_string(Kmer* kmer) {
	char str[32];
	//! elog(INFO, "kmer value (kmer_value_to_string): %lu", kmer -> value);
	for (uint8_t i = 0; i < kmer -> k; i++) {
        uint8_t shift = (kmer -> k - i - 1) * 2;
        uint8_t nucleotide = (kmer -> value >> shift) & 0b11;

		str[i] = BINARY_TO_NUCLEOTIDE[nucleotide];
    }
	str[kmer->k] = '\0';
	return psprintf("%s", str);
}
/**
 * @brief Checks if a K-mer starts with a prefix.
 * 
 * @param kmer The K-mer to check.
 * @param prefix The prefix to check.
 * @return True if the K-mer starts with the prefix, false otherwise.
 */
bool internal_kmer_startswith(Kmer* kmer, Kmer* prefix) {
	if (kmer -> k < prefix -> k) {
		return false;
	}
	uint64_t extracted_from_kmer = kmer -> value >> (kmer -> k - prefix -> k) * 2;
	return extracted_from_kmer == prefix -> value;
}

/**
 * @brief Get the first k nucleotides of a K-mer.
 * 
 * @param kmer The K-mer.
 * @param k The number of nucleotides to get.
 * @return The K-mer with the first k nucleotides.
 */
Kmer* get_first_k_nucleotides(Kmer* kmer, uint8_t k) {
    if (k > kmer->k) {
        ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE), errmsg("k cannot be greater than the K-mer size")));
    }
    Kmer* first_kmer = palloc0(sizeof(Kmer));
    first_kmer->k = k;
    first_kmer->value = kmer->value >> (2 * (kmer->k - k));
    return first_kmer;
}

/**
 * @brief Get the last k nucleotides of a K-mer.
 * 
 * @param kmer The K-mer.
 * @param k The number of nucleotides to get.
 * @return The K-mer with the last k nucleotides.
 */
Kmer* get_last_k_nucleotides(Kmer* kmer, uint8_t k) {
    Kmer* last_kmer = palloc0(sizeof(Kmer));
    last_kmer->k = k;
    last_kmer->value = kmer->value & ((1 << (2 * k)) - 1);
    return last_kmer;
}

/**
 * @brief Function to get the common prefix length of 2 K-mers.
 * 
 * @param kmer1 The first K-mer.
 * @param kmer2 The second K-mer.
 * @return The common prefix length.
 */
uint8_t get_common_prefix_len(Kmer* kmer1, Kmer* kmer2) {
    uint64_t kmer1_value = kmer1->value; // Copy the value to avoid modifying the original
    uint64_t kmer2_value = kmer2->value; // Copy the value to avoid modifying the original

    int size_diff = kmer1->k - kmer2->k;
    uint8_t prefix_len;

    /*
     * If the K-mers have different sizes, we need to align them by shifting the larger one.
     * We then compare the values to find the common prefix length.
     */
    if (size_diff > 0) {
        prefix_len = kmer2->k;
        kmer1_value >>= 2 * size_diff;
    } else if (size_diff < 0) {
        prefix_len = kmer1->k;
        kmer2_value >>= 2 * -size_diff;
    } else {
        prefix_len = kmer1->k;
    }
    uint64_t xor = kmer1_value ^ kmer2_value;
    while (xor != 0) {                            // if prefix_len is 0, xor will be 0
        xor >>= 2;
        prefix_len--;
    }
    return prefix_len;
}

/**
 * @brief Compare the n first nucleotides of two K-mers.
 * 
 * @param kmer1 The first K-mer.
 * @param kmer2 The second K-mer.
 * @param n The number of nucleotides to compare.
 * @return The comparison result, -1 if kmer1 < kmer2, 0 if kmer1 == kmer2, 1 if kmer1 > kmer2.
 */
int compare_kmers(Kmer* kmer1, Kmer* kmer2, uint8_t n) {
    Kmer* first_kmer1 = get_first_k_nucleotides(kmer1, n);
    Kmer* first_kmer2 = get_first_k_nucleotides(kmer2, n);

    int result = 0;
    if (first_kmer1->value < first_kmer2->value) {
        result = -1;
    } else if (first_kmer1->value > first_kmer2->value) {
        result = 1;
    }
    pfree(first_kmer1);
    pfree(first_kmer2);
    return result;
}


/* ************************************************************************** */

/**
 * @brief Postgres input function for K-mer.
 * 
 * @param str The input string.
 * @return The K-mer object created from the input string.
 */
PG_FUNCTION_INFO_V1(kmer_in);
Datum kmer_in(PG_FUNCTION_ARGS) {
	char *str = PG_GETARG_CSTRING(0);
	Kmer* kmer = kmer_parse(str);
	PG_RETURN_KMER_P(kmer);
}

/**
 * @brief Postgres output function for K-mer.
 * 
 * @param kmer The K-mer object.
 * @return The string representation of the K-mer.
 */
PG_FUNCTION_INFO_V1(kmer_out);
Datum kmer_out(PG_FUNCTION_ARGS) {
	Kmer *kmer = PG_GETARG_KMER_P(0);
	char *str = kmer_value_to_string(kmer);
	PG_FREE_IF_COPY(kmer, 0);
	PG_RETURN_CSTRING(str);
}

/**
 * @brief Postgres receive function for K-mer.
 * 
 * @param buf The bytea representation of the K-mer.
 * @return The K-mer object created from the bytea.
 */
PG_FUNCTION_INFO_V1(kmer_recv);
Datum kmer_recv(PG_FUNCTION_ARGS) {
	StringInfo buf = (StringInfo) PG_GETARG_POINTER(0);
	Kmer *kmer = palloc0(sizeof(Kmer));
	kmer -> value = pq_getmsgint64(buf);
	kmer -> k = pq_getmsgint(buf, 8);
	PG_RETURN_KMER_P(kmer);
}

/**
 * @brief Postgres send function for K-mer.
 * 
 * @param kmer The K-mer object.
 * @return The bytea representation of the K-mer.
 */
PG_FUNCTION_INFO_V1(kmer_send);
Datum kmer_send(PG_FUNCTION_ARGS) {
	Kmer *kmer = PG_GETARG_KMER_P(0);
	StringInfoData buf;
	pq_begintypsend(&buf);
	pq_sendint64(&buf, kmer -> value);
	pq_sendint8(&buf, kmer -> k);
	PG_FREE_IF_COPY(kmer, 0);
	PG_RETURN_BYTEA_P(pq_endtypsend(&buf));
}

/**
 * @brief Postgres cast function from text to K-mer.
 * 
 * @param txt The text to cast.
 * @return The K-mer object created from the text.
 */
PG_FUNCTION_INFO_V1(kmer_cast_from_text);
Datum kmer_cast_from_text(PG_FUNCTION_ARGS) {
	text *txt = PG_GETARG_TEXT_P(0);
	char *str = DatumGetCString(DirectFunctionCall1(textout,
	             PointerGetDatum(txt)));
	Kmer* kmer = kmer_parse(str);
	PG_FREE_IF_COPY(txt, 0);
	PG_RETURN_KMER_P(kmer);
}

/**
 * @brief Postgres cast function from K-mer to text.
 * 
 * @param kmer The K-mer to cast.
 * @return The text representation of the K-mer.
 */
PG_FUNCTION_INFO_V1(kmer_cast_to_text);
Datum kmer_cast_to_text(PG_FUNCTION_ARGS) {
	Kmer* kmer  = PG_GETARG_KMER_P(0);
	text* out = (text *)DirectFunctionCall1(textin, 
						PointerGetDatum(kmer_value_to_string(kmer)));
	PG_FREE_IF_COPY(kmer, 0);
	PG_RETURN_TEXT_P(out);
}

/**
 * @brief Postgres length function for K-mer.
 * 
 * @param kmer The K-mer to get the length of.
 * @return The length of the K-mer.
 */
PG_FUNCTION_INFO_V1(kmer_length);
Datum kmer_length(PG_FUNCTION_ARGS) {
	Kmer* kmer  = PG_GETARG_KMER_P(0);
	uint8_t length = kmer -> k;
	PG_FREE_IF_COPY(kmer, 0);
	PG_RETURN_CHAR(length);
}

/**
 * @brief Postgres function to check if two K-mers are equal.
 * 
 * @param a The first K-mer.
 * @param b The second K-mer.
 * @return True if the K-mers are equal, false otherwise.
 */
PG_FUNCTION_INFO_V1(kmer_eq);
Datum kmer_eq(PG_FUNCTION_ARGS) {
	Kmer* a = PG_GETARG_KMER_P(0);
	Kmer* b = PG_GETARG_KMER_P(1);
	bool result = KMER_EQUAL(a, b);
	PG_FREE_IF_COPY(a, 0);
	PG_FREE_IF_COPY(b, 1);
	PG_RETURN_BOOL(result);
}

/**
 * @brief Postgres function to check if a K-mer starts with a prefix.
 * 
 * @param kmer The K-mer to check.
 * @param prefix The prefix to check.
 * @return True if the K-mer starts with the prefix, false otherwise.
 */
PG_FUNCTION_INFO_V1(kmer_startswith);
Datum kmer_startswith(PG_FUNCTION_ARGS) {
	Kmer *kmer = PG_GETARG_KMER_P(0);
	Kmer *prefix = PG_GETARG_KMER_P(1);
	bool result = internal_kmer_startswith(kmer, prefix);
	PG_FREE_IF_COPY(kmer, 0);
	PG_FREE_IF_COPY(prefix, 1);
	PG_RETURN_BOOL(result);
}

/* Kmer Hash operators */

/**
 * @brief Postgres function to get the hash value of a K-mer.
 * 
 * @param kmer The K-mer to get the hash value of.
 * @return The hash value of the K-mer.
 */
PG_FUNCTION_INFO_V1(kmer_hash);
Datum kmer_hash(PG_FUNCTION_ARGS) {
	Kmer* kmer = PG_GETARG_KMER_P(0);

	uint64_t kmer_value = kmer->value << (kmer->k * 2);
	uint64_t removed_bits = kmer->value >> (64 - kmer->k * 2);
	uint64_t hash_input = kmer_value ^ removed_bits;

	int32 hash = hash_any((unsigned char *) &hash_input, sizeof(hash_input));

    PG_RETURN_INT32(hash);
}