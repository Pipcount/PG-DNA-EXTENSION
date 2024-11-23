#ifndef KMEA_H  // include guard
#define KMEA_H

#include "postgres.h"
#include "fmgr.h"

#include "libpq/pqformat.h"
#include "utils/fmgrprotos.h"
#include "utils/elog.h"
#include "utils/varlena.h"
#include "varatt.h"

// Define macros for Kmer because we use a struct to represent a Kmer
#define DatumGetKmerP(X)  ((Kmer *) DatumGetPointer(X))
#define KmerPGetDatum(X)  PointerGetDatum(X)
#define PG_GETARG_KMER_P(n) DatumGetKmerP(PG_GETARG_DATUM(n))
#define PG_RETURN_KMER_P(x) return KmerPGetDatum(x)

// Define macros for Qkmer because we use a struct to represent a Qkmer
#define DatumGetQkmerP(X)  ((Qkmer *) DatumGetPointer(X))
#define QkmerPGetDatum(X)  PointerGetDatum(X)
#define PG_GETARG_QKMER_P(n) DatumGetQkmerP(PG_GETARG_DATUM(n))
#define PG_RETURN_QKMER_P(x) return QkmerPGetDatum(x)

/** 
 * @typedef DNA
 * @brief Type used to store a DNA sequence.
 */
typedef bytea DNA;

/**
 * @typedef Kmer
 * @brief Structure used to store a K-mer.
 */
typedef struct Kmer {
	uint64_t value;	  /**< The value of the K-mer */
	uint8_t k;        /**< The length of the K-mer */
} Kmer;

/**
 * @typedef Qkmer
 * @brief Structure used to store a Q-kmer.
 */
typedef struct Qkmer
{
    uint64_t ac;    /**< The value of the A/C part of the Q-kmer */
    uint64_t gt;    /**< The value of the G/T part of the Q-kmer */
    uint8_t  k;     /**< The length of the Q-kmer */
} Qkmer;

/**
 * @brief LUT to convert a nucleotide to a 2-bit representation.
 * 
 * The LUT is indexed by the ASCII value of the nucleotide minus the ASCII value of 'A'.
 * Encoding for non-nucleotide characters is 0. This overlaps with the encoding for 'A' but should not be a problem if we check the input.
 */
static const uint8_t NUCLEOTIDE_TO_BINARY[26] = {
    0b00,  // A
    0,  // B
    0b01,  // C
    0,  // D
    0,  // E
    0,  // F
    0b10,  // G
    0,  // H
    0,  // I
    0,  // J
    0,  // K
    0,  // L
    0,  // M
    0,  // N
    0,  // O
    0,  // P
    0,  // Q
    0,  // R
    0,  // S
    0b11,  // T
    0,  // U
    0,  // V
    0,  // W
    0,  // X
    0,  // Y
    0   // Z
};

/**
 * @brief LUT to convert a 2-bit representation to a nucleotide.
 */
static const char BINARY_TO_NUCLEOTIDE[4] = {
    'A',  // 00
    'C',  // 01
    'G',  // 10
    'T'   // 11
};

/**
 * @brief Encodes a nucleotide to a 2-bit representation and adds it to a binary.
 * 
 * @param binary An unsigned integer to which the nucleotide will be added.
 * @param nucleotide The nucleotide to encode and add.
 * @return void
 */
#define add_nucleotide_to_uint(binary, nucleotide) \
    do { \
        if (!strchr("ACGT", toupper(nucleotide))) { \
            ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION), \
                errmsg("invalid nucleotide"))); \
        } \
        binary = (binary << 2) | NUCLEOTIDE_TO_BINARY[toupper(nucleotide) - 'A']; \
    } while (0)

#endif