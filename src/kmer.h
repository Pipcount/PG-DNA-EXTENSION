#ifndef KMER_H
#define KMER_H

#include "kmea.h"
#include <stdint.h>
#include <stdio.h>
#include <string.h>

static Kmer* make_kmer(const char *str, uint8_t length);
static Kmer* kmer_parse(const char* str);
static char* kmer_value_to_string(Kmer* kmer);

#endif