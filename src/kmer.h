#ifndef KMER_H
#define KMER_H

#include "kmea.h"
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "access/hash.h"

// Macro to check if two Kmers are equal
#define KMER_EQUAL(kmer1, kmer2) ((kmer1 -> value == kmer2 -> value) && (kmer1 -> k == kmer2 -> k))

bool internal_kmer_startswith(Kmer* kmer, Kmer* prefix);
char* kmer_value_to_string(Kmer* kmer);

Kmer* get_first_k_nucleotides(Kmer* kmer, uint8_t k);
Kmer* get_last_k_nucleotides(Kmer* kmer, uint8_t k);

uint8_t get_common_prefix_len(Kmer* kmer1, Kmer* kmer2);
int compare_kmers(Kmer* kmer1, Kmer* kmer2, uint8_t n);

#endif