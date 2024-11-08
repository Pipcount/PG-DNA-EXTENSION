#ifndef DNA_H
#define DNA_H

#include "kmea.h"
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "funcapi.h"

static DNA* make_dna(const char* str, uint32_t length);
static DNA* dna_parse(const char* str);
static char* dna_to_string(DNA* dna);

static void store_last_byte_length(DNA* dna, uint32_t length);
static uint32_t get_dna_sequence_length(DNA* dna);

#endif