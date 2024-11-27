#ifndef KMER_H
#define KMER_H

#include "kmea.h"
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "access/hash.h"

// Macro to check if two Kmers are equal
#define KMER_EQUAL(kmer1, kmer2) ((kmer1 -> value == kmer2 -> value) && (kmer1 -> k == kmer2 -> k))

#endif