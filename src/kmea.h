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


#endif