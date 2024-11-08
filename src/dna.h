#ifndef DNA_H
#define DNA_H

#include "kmea.h"
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "funcapi.h"

typedef struct KmerGeneratorState {
    uint8_t kmer_length;
    uint32_t length;
    uint8_t *byte_ptr;
    uint8_t nucleotide_ctr;
    uint32_t current_pos;
} KmerGeneratorState;

static DNA* make_dna(const char* str, uint32_t length);
static DNA* dna_parse(const char* str);
static char* dna_to_string(DNA* dna);

static void store_last_byte_length(DNA* dna, uint32_t length);
static uint32_t get_dna_sequence_length(DNA* dna);
void init_kmer_generator_state(KmerGeneratorState *state, DNA *dna, FuncCallContext *funcctx);

#endif

