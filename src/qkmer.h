#ifndef QKMER_H
#define QKMER_H

#include "kmea.h"

Qkmer* get_first_k_nucleotides_qkmer(Qkmer* qkmer, uint8_t k);

bool qkmer_contains_internal(Qkmer* qkmer, Kmer* kmer);

#endif