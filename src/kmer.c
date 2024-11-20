#include "kmer.h"

PG_MODULE_MAGIC;

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
		switch (toupper(c)) {
			case 'A': kmer -> value = (kmer -> value << 2) | 0b00; break; // TODO: Use LUT for this
			case 'C': kmer -> value = (kmer -> value << 2) | 0b01; break;
			case 'G': kmer -> value = (kmer -> value << 2) | 0b10; break;
			case 'T': kmer -> value = (kmer -> value << 2) | 0b11; break;
			default:
				ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
				errmsg("invalid nucleotide")));
				break;
		}
	}
	//! elog(INFO, "kmer value after make_kmer: %lu", kmer -> value);
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
static char* kmer_value_to_string(Kmer* kmer) {
	char str[32];
	//! elog(INFO, "kmer value (kmer_value_to_string): %lu", kmer -> value);
	for (uint8_t i = 0; i < kmer -> k; i++) {
        uint8_t shift = (kmer -> k - i - 1) * 2;
        uint8_t nucleotide = (kmer -> value >> shift) & 0b11;

        if (nucleotide == 0b00) {       // TODO: Use LUT for this
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
/**
 * @brief Checks if a K-mer starts with a prefix.
 * 
 * @param kmer The K-mer to check.
 * @param prefix The prefix to check.
 * @return True if the K-mer starts with the prefix, false otherwise.
 */
static bool internal_kmer_startswith(Kmer* kmer, Kmer* prefix) {
	if (kmer -> k < prefix -> k) {
		ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
	  	errmsg("kmer should not be shorter than prefix")));
		return false;
	}
	uint64_t extracted_from_kmer = kmer -> value >> (kmer -> k - prefix -> k) * 2;
	return extracted_from_kmer == prefix -> value;
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
	//! elog(INFO, "kmer value (kmer_in): %lu", kmer -> value);
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
	//! elog(INFO, "kmer value (kmer_out): %lu, size %d", kmer -> value, kmer -> k);
	char *str = kmer_value_to_string(kmer);
	//! elog(INFO, "kmer string (kmer_out): %sblablabla", str);
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
	//! elog(INFO, "kmer_recv");
	StringInfo buf = (StringInfo) PG_GETARG_POINTER(0);
	Kmer *kmer = palloc(sizeof(Kmer));
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
	//! elog(INFO, "kmer_send");
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
Datum
kmer_cast_from_text(PG_FUNCTION_ARGS)
{
  text *txt = PG_GETARG_TEXT_P(0);
  char *str = DatumGetCString(DirectFunctionCall1(textout,
               PointerGetDatum(txt)));
  PG_RETURN_KMER_P(kmer_parse(str));
}

/**
 * @brief Postgres cast function from K-mer to text.
 * 
 * @param kmer The K-mer to cast.
 * @return The text representation of the K-mer.
 */
PG_FUNCTION_INFO_V1(kmer_cast_to_text);
Datum
kmer_cast_to_text(PG_FUNCTION_ARGS)
{
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
Datum
kmer_length(PG_FUNCTION_ARGS)
{
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

