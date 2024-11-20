#include "qkmer.h"

PG_MODULE_MAGIC;


/** 
 * @brief Creates a Q-kmer from a string.
 * 
 * @param str The string representing the Q-kmer.
 * @param length The length of the Q-kmer.
 * @return A pointer to the created Q-kmer.
 */
static Qkmer* make_qkmer_from_str(const char* str, uint8_t length) {
    Qkmer* qkmer = palloc0(sizeof(Qkmer));
    qkmer -> k = length;
    qkmer -> ac = 0;
    qkmer -> gt = 0;

    for (uint8_t i = 0; i < length; i++) {
        char c = str[i];
        switch (toupper(c)) {
            case 'A': qkmer -> ac = (qkmer -> ac << 2) | 0b10; qkmer -> gt = (qkmer -> gt << 2) | 0b00; break; // TODO: Use LUT for this
            case 'C': qkmer -> ac = (qkmer -> ac << 2) | 0b01; qkmer -> gt = (qkmer -> gt << 2) | 0b00; break;
            case 'G': qkmer -> ac = (qkmer -> ac << 2) | 0b00; qkmer -> gt = (qkmer -> gt << 2) | 0b10; break;
            case 'T': qkmer -> ac = (qkmer -> ac << 2) | 0b00; qkmer -> gt = (qkmer -> gt << 2) | 0b01; break;
            case 'W': qkmer -> ac = (qkmer -> ac << 2) | 0b10; qkmer -> gt = (qkmer -> gt << 2) | 0b01; break;
            case 'S': qkmer -> ac = (qkmer -> ac << 2) | 0b01; qkmer -> gt = (qkmer -> gt << 2) | 0b10; break;
            case 'M': qkmer -> ac = (qkmer -> ac << 2) | 0b11; qkmer -> gt = (qkmer -> gt << 2) | 0b00; break;
            case 'K': qkmer -> ac = (qkmer -> ac << 2) | 0b00; qkmer -> gt = (qkmer -> gt << 2) | 0b11; break;
            case 'R': qkmer -> ac = (qkmer -> ac << 2) | 0b10; qkmer -> gt = (qkmer -> gt << 2) | 0b10; break;
            case 'Y': qkmer -> ac = (qkmer -> ac << 2) | 0b01; qkmer -> gt = (qkmer -> gt << 2) | 0b01; break;
            case 'B': qkmer -> ac = (qkmer -> ac << 2) | 0b01; qkmer -> gt = (qkmer -> gt << 2) | 0b11; break;
            case 'D': qkmer -> ac = (qkmer -> ac << 2) | 0b10; qkmer -> gt = (qkmer -> gt << 2) | 0b11; break;
            case 'H': qkmer -> ac = (qkmer -> ac << 2) | 0b11; qkmer -> gt = (qkmer -> gt << 2) | 0b01; break;
            case 'V': qkmer -> ac = (qkmer -> ac << 2) | 0b11; qkmer -> gt = (qkmer -> gt << 2) | 0b10; break;
            case 'N': qkmer -> ac = (qkmer -> ac << 2) | 0b11; qkmer -> gt = (qkmer -> gt << 2) | 0b11; break;
            default:
                ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                errmsg("invalid nucleotide")));
                break;
        }
    }
    return qkmer;
}

/**
 * @brief Creates a Q-kmer from a K-mer.
 * 
 * @param kmer The K-mer to create the Q-kmer from.
 * @return A pointer to the created Q-kmer.
 */
static Qkmer* make_qkmer_from_kmer(Kmer* kmer) {
    return NULL; // TODO
}
/**
 * @brief Parses a Q-kmer from a string.
 * 
 * @param str The string representing the Q-kmer to parse.
 * @return A pointer to the Q-kmer created from the string.
 */
static Qkmer* qkmer_parse(const char* str) {
    uint8_t length = strlen(str);

	if (length > 32) {
		ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
      	errmsg("qkmer should not exceed 32 nucleotides")));
	} else if (length == 0) {
		ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
	  	errmsg("qkmer should not be empty")));
	}
	return make_qkmer_from_str(str, length);
}
/* ************************************************************************** */

/**
 * @brief Postgres input function for Q-kmer.
 * 
 * @param str The input string.
 * @return The Q-kmer object created from the input string.
 */
PG_FUNCTION_INFO_V1(qkmer_in);
Datum qkmer_in(PG_FUNCTION_ARGS) {
    char* str = PG_GETARG_CSTRING(0);
    Qkmer* qkmer = qkmer_parse(str);
    PG_RETURN_QKMER_P(qkmer);
}

/**
 * @brief Postgres output function for Q-kmer.
 * 
 * @param qkmer The Q-kmer object.
 * @return The string representation of the Q-kmer.
 */
PG_FUNCTION_INFO_V1(qkmer_out);
Datum qkmer_out(PG_FUNCTION_ARGS) {

}

/**
 * @brief Postgres receive function for Q-kmer.
 * 
 * @param buf The bytea representation of the Q-kmer.
 * @return The Q-kmer object created from the bytea.
 */
PG_FUNCTION_INFO_V1(qkmer_send);
Datum qkmer_send(PG_FUNCTION_ARGS) {

}

/**
 * @brief Postgres receive function for Q-kmer.
 * 
 * @param buf The bytea representation of the Q-kmer.
 * @return The Q-kmer object created from the bytea.
 */
PG_FUNCTION_INFO_V1(qkmer_recv);
Datum qkmer_recv(PG_FUNCTION_ARGS) {

}