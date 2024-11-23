#ifndef KMEA_H  // include guard
#define KMEA_H

#include "postgres.h"
#include "fmgr.h"

#include "libpq/pqformat.h"
#include "utils/fmgrprotos.h"
#include "utils/elog.h"
#include "utils/varlena.h"
#include "varatt.h"

/** 
 * @typedef DNA
 * @brief Type used to store a DNA sequence.
 */
typedef struct DNA {
    uint8_t last_byte_length;       /**< The length of the last byte of the DNA sequence */
    bytea sequence;                 /**< The byte array representing the DNA sequence */
} DNA;

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
    uint32_t k;     /**< The length of the Q-kmer */
} Qkmer;

// Define macros for DNA because we use a struct to represent a DNA sequence
#define DatumGetDNAP(X)  ((DNA *) DatumGetPointer(X))
#define DNAPGetDatum(X)  PointerGetDatum(X)
#define PG_GETARG_DNA_P(n) DatumGetDNAP(PG_GETARG_DATUM(n))
#define PG_RETURN_DNA_P(x) return DNAPGetDatum(x)

// Define macros for Kmer because we use a struct to represent a Kmer
#define DatumGetKmerP(X)  ((Kmer *) DatumGetPointer(X))
#define KmerPGetDatum(X)  PointerGetDatum(X)
#define PG_GETARG_KMER_P(n) DatumGetKmerP(PG_GETARG_DATUM(n))
#define PG_RETURN_KMER_P(x) return KmerPGetDatum(x)

// Define macros for Qkmer because we use a struct to represent a Qkmer
#define DatumGetQkmerP(X)  ((Qkmer *) DatumGetPointer(X))
#define QkmerPGetDatum(X)  PointerGetDatum(X)
#define PG_GETARG_QKMER_P(n) DatumGetGkmerP(PG_GETARG_DATUM(n))
#define PG_RETURN_QKMER_P(x) return QkmerPGetDatum(x)

#endif