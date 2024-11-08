#ifndef KMEA_H  // include guard
#define KMEA_H

#include "postgres.h"
#include "fmgr.h"

#include "libpq/pqformat.h"
#include "utils/fmgrprotos.h"
#include "utils/elog.h"
#include "utils/varlena.h"
#include "varatt.h"

typedef bytea DNA;

typedef struct Kmer {
	uint64_t value;	// kmer value
	uint8_t k;  // kmer length
} Kmer;

// Define macros for Kmer because we use a struct to represent a Kmer
#define DatumGetKmerP(X)  ((Kmer *) DatumGetPointer(X))
#define KmerPGetDatum(X)  PointerGetDatum(X)
#define PG_GETARG_KMER_P(n) DatumGetKmerP(PG_GETARG_DATUM(n))
#define PG_RETURN_KMER_P(x) return KmerPGetDatum(x)


#endif