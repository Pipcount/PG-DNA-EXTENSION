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
 * @brief Converts a Q-kmer to a string.
 * 
 * @param qkmer The Q-kmer to convert.
 * @return The string representation of the Q-kmer.
 */
static char* qkmer_value_to_string(Qkmer* qkmer) {
    char str[32];
    for (uint8_t i = 0; i < qkmer -> k; i++) {
        uint8_t shift = (qkmer -> k - i - 1) * 2;
        uint8_t ac_nucleotide = (qkmer -> ac >> shift) & 0b11;
        uint8_t gt_nucleotide = (qkmer -> gt >> shift) & 0b11;

        uint8_t nucleotide = (ac_nucleotide << 2) | gt_nucleotide;
        switch (nucleotide) {
            case 0b1000: str[i] = 'A'; break;        // TODO: Use LUT for this
            case 0b0100: str[i] = 'C'; break;
            case 0b0010: str[i] = 'G'; break;
            case 0b0001: str[i] = 'T'; break;
            case 0b1001: str[i] = 'W'; break;
            case 0b0110: str[i] = 'S'; break;
            case 0b1100: str[i] = 'M'; break;
            case 0b0011: str[i] = 'K'; break;
            case 0b1010: str[i] = 'R'; break;
            case 0b0101: str[i] = 'Y'; break;
            case 0b0111: str[i] = 'B'; break;
            case 0b1011: str[i] = 'D'; break;
            case 0b1101: str[i] = 'H'; break;
            case 0b1110: str[i] = 'V'; break;
            case 0b1111: str[i] = 'N'; break;
            default:
                ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                                errmsg("invalid binary value")));
        }
    }
    str[qkmer -> k] = '\0';
    return psprintf("%s", str);
}

/**
 * @brief Creates a Q-kmer from a K-mer.
 * 
 * @param kmer The K-mer to create the Q-kmer from.
 * @return A pointer to the created Q-kmer.
 */
static Qkmer* make_qkmer_from_kmer(Kmer* kmer) {
    Qkmer* qkmer = palloc0(sizeof(Qkmer));
    qkmer -> k = kmer -> k;

    for (uint8_t i = 0; i < kmer -> k; i++) {
        uint8_t shift = (kmer -> k - i - 1) * 2;
        uint8_t nucleotide = (kmer -> value >> shift) & 0b11;

        switch (nucleotide) {       // TODO: Use LUT for this
            case 0b00:
            qkmer -> ac = (qkmer -> ac << 2) | 0b10;
            qkmer -> gt = (qkmer -> gt << 2) | 0b00;
            break;
            case 0b01:
            qkmer -> ac = (qkmer -> ac << 2) | 0b01;
            qkmer -> gt = (qkmer -> gt << 2) | 0b00;
            break;
            case 0b10:
            qkmer -> ac = (qkmer -> ac << 2) | 0b00;
            qkmer -> gt = (qkmer -> gt << 2) | 0b10;
            break;
            case 0b11:
            qkmer -> ac = (qkmer -> ac << 2) | 0b00;
            qkmer -> gt = (qkmer -> gt << 2) | 0b01;
            break;
            default:
            ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                    errmsg("invalid nucleotide for kmer")));
            break;
        }
    }
    return qkmer;
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
    Qkmer* qkmer = PG_GETARG_QKMER_P(0);
    char* str = qkmer_value_to_string(qkmer);
    PG_FREE_IF_COPY(qkmer, 0);
    PG_RETURN_CSTRING(str);
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

/**
 * @brief Postgres cast function from text to Q-kmer.
 * 
 * @param txt The text to cast.
 * @return The Q-kmer object created from the text.
 */
PG_FUNCTION_INFO_V1(qkmer_cast_from_text);
Datum qkmer_cast_from_text(PG_FUNCTION_ARGS) {
    text* txt = PG_GETARG_TEXT_P(0);
    char* str = DatumGetCString(DirectFunctionCall1(textout, PointerGetDatum(txt)));
    Qkmer* qkmer = qkmer_parse(str);
    PG_FREE_IF_COPY(txt, 0);
    PG_RETURN_QKMER_P(qkmer);
}

/**
 * @brief Postgres cast function from Q-kmer to text.
 * 
 * @param qkmer The Q-kmer to cast.
 * @return The text representation of the Q-kmer.
 */
PG_FUNCTION_INFO_V1(qkmer_cast_to_text);
Datum qkmer_cast_to_text(PG_FUNCTION_ARGS) {
    Qkmer* qkmer = PG_GETARG_QKMER_P(0);
    text* out = (text *)DirectFunctionCall1(textin, PointerGetDatum(qkmer_value_to_string(qkmer)));
    PG_FREE_IF_COPY(qkmer, 0);
    PG_RETURN_TEXT_P(out);
}

/**
 * @brief Checks if a Q-kmer matches a kmer
 * 
 * @param qkmer The Q-kmer to check.
 * @param kmer The K-mer to check.
 * @return True if the Q-kmer matches the K-mer, false otherwise.
 */
PG_FUNCTION_INFO_V1(qkmer_contains);
Datum qkmer_contains(PG_FUNCTION_ARGS) {
    Qkmer* qkmer = PG_GETARG_QKMER_P(0);
    Kmer* kmer = PG_GETARG_KMER_P(1);
    
    if (qkmer -> k != kmer -> k) {
        PG_RETURN_BOOL(false);
    }

    Qkmer* qkmer_from_kmer = make_qkmer_from_kmer(kmer);
    bool result = (qkmer_from_kmer -> ac & qkmer -> ac) == qkmer_from_kmer -> ac &&
                  (qkmer_from_kmer -> gt & qkmer -> gt) == qkmer_from_kmer -> gt;
    PG_FREE_IF_COPY(qkmer, 0);
    PG_FREE_IF_COPY(kmer, 1);
    PG_RETURN_BOOL(result);
}

/**
 * @brief Returns the length of a Q-kmer.
 * 
 * @param qkmer The Q-kmer to get the length of.
 * @return The length of the Q-kmer.
 */
PG_FUNCTION_INFO_V1(qkmer_length);
Datum qkmer_length(PG_FUNCTION_ARGS) {
    Qkmer* qkmer = PG_GETARG_QKMER_P(0);
    uint8_t length = qkmer -> k;
    PG_FREE_IF_COPY(qkmer, 0);
    PG_RETURN_CHAR(length);
}