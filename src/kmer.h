#ifndef KMER_H
#define KMER_H

#include "kmea.h"
#include <stdint.h>
#include <stdio.h>
#include <string.h>


// Define macros for Kmer because we use a struct to represent a Kmer
#define DatumGetKmerP(X)  ((Kmer *) DatumGetPointer(X))
#define KmerPGetDatum(X)  PointerGetDatum(X)
#define PG_GETARG_KMER_P(n) DatumGetKmerP(PG_GETARG_DATUM(n))
#define PG_RETURN_KMER_P(x) return KmerPGetDatum(x)

static Kmer* make_kmer(const char *str, uint8_t length);
static Kmer* kmer_parse(const char* str);
static char* kmer_value_to_string(Kmer* kmer);

#endif