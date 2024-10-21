#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "postgres.h"

PG_MODULE_MAGIC;

typedef int64_t kmer_t;

static inline kmer_t kmerFromStr(const char* str) {

	int str_len = strlen(str);
	if (str_len > 32) {
		perror("kmerFromStr: kmer is too long");
	}

	kmer_t kmer = 0;
	

}
