#ifndef KMER_H
#define KMER_H

#include "kmea.h"
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define KMER_EQUAL(kmer1, kmer2) ((kmer1 -> value == kmer2 -> value) && (kmer1 -> k == kmer2 -> k))

static Kmer* make_kmer(const char *str, uint8_t length);
static Kmer* kmer_parse(const char* str);
static char* kmer_value_to_string(Kmer* kmer);
static bool internal_kmer_startswith(Kmer* a, Kmer* b);

#endif