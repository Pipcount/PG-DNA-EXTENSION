#include "qkmer.h"

/**
 * @brief LUT to convert a IUPAC nucleotide code to a 4-bit representation.
 * 
 * The LUT is indexed by the ASCII value of the nucleotide minus the ASCII value of 'A'.
 * Encoding for non IUPAC nucleotide characters is {0b00, 0b00}. This needs to be checked when using the LUT.
 */
static const uint8_t IUPAC_CODE_TO_BINARY[26][2] = {
    {0b10, 0b00},  // A
    {0b01, 0b11},  // B
    {0b01, 0b00},  // C
    {0b10, 0b11},  // D
    {0b00, 0b00},  // E
    {0b00, 0b00},  // F
    {0b00, 0b10},  // G
    {0b11, 0b01},  // H
    {0b00, 0b00},  // I
    {0b00, 0b00},  // J
    {0b00, 0b11},  // K
    {0b00, 0b00},  // L
    {0b11, 0b00},  // M
    {0b11, 0b11},  // N
    {0b00, 0b00},  // O
    {0b00, 0b00},  // P
    {0b00, 0b00},  // Q
    {0b10, 0b10},  // R
    {0b01, 0b10},  // S
    {0b00, 0b01},  // T
    {0b00, 0b00},  // U
    {0b11, 0b10},  // V
    {0b10, 0b01},  // W
    {0b00, 0b00},  // X
    {0b01, 0b01},  // Y
    {0b00, 0b00}   // Z     
};

/**
 * @brief LUT to convert a 4-bit representation of an IUPAC nucleotide code to a character representation.
 * 
 * If the 4-bit representation is 0, the nucleotide is invalid and the character is '$'.
 * This should never happen but is handled here just in case.
 */
static const char BINARY_TO_IUPAC_CODE[16] = {
    '$',  // 0000
    'T',   // 0001
    'G',   // 0010
    'K',   // 0011
    'C',   // 0100
    'Y',   // 0101
    'S',   // 0110
    'B',   // 0111
    'A',   // 1000
    'W',   // 1001
    'R',   // 1010
    'D',   // 1011
    'M',   // 1100
    'H',   // 1101
    'V',   // 1110
    'N'    // 1111
};

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
        uint8_t index = toupper(c) - 'A';
        if (index < 0 || index > 25 || (IUPAC_CODE_TO_BINARY[index][0] == 0b00 && IUPAC_CODE_TO_BINARY[index][1] == 0b00)) {
            ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                            errmsg("invalid nucleotide")));
        }
        qkmer -> ac = (qkmer -> ac << 2) | IUPAC_CODE_TO_BINARY[index][0];
        qkmer -> gt = (qkmer -> gt << 2) | IUPAC_CODE_TO_BINARY[index][1];

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
        char c = BINARY_TO_IUPAC_CODE[nucleotide];
        if (c == '$') {                                             // This is maybe not the best way to handle this
            ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION), // This should also never happen
                            errmsg("invalid nucleotide")));
        }
        str[i] = c;
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

    // Masks for extracting 2-bit pairs
    const uint64_t odd_position_pair_mask = 0x3333333333333333; // Binary: 00110011... (two 1s in each pair)
    const uint64_t even_position_pair_mask = 0xCCCCCCCCCCCCCCCC; // Binary: 11001100... (two 0s in each pair)
    const uint64_t zero_one_mask = 0x5555555555555555; // Binary: 01010101... 
    const uint64_t one_zero_mask = 0xAAAAAAAAAAAAAAAA; // Binary: 10101010... 

    // Isolate T (11):
    uint64_t t_odd = kmer->value & odd_position_pair_mask;
    uint64_t t_even = kmer->value & even_position_pair_mask;
    uint64_t t = t_odd & (t_odd >> 1) | t_even & (t_even >> 1);

    // Isolate A (00): 
    uint64_t a_odd = ~kmer->value & odd_position_pair_mask; 
    uint64_t a_even = ~kmer->value & even_position_pair_mask;   
    uint64_t a = a_odd & (a_odd >> 1) | a_even & (a_even >> 1);
    a = a << 1;                                                     // to have 10 instead of 01

    // Isolate C (01):
    uint64_t c = kmer->value & zero_one_mask & ~t;

    // Isolate G (10):
    uint64_t g = (kmer->value & one_zero_mask) >> 1 & ~t;
    g = g << 1;                                                     // to have 10 instead of 01

    uint64_t length_mask = (1ULL << (kmer -> k * 2)) - 1;

    qkmer -> ac = (a | c) & length_mask; // Mask to get rid of the bits that are not part of the k-mer (is necessary here)
    qkmer -> gt = (g | t) & length_mask; // Mask to get rid of the bits that are not part of the k-mer (should not shange anything here)
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