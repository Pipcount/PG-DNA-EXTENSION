#ifndef QKMER_H
#define QKMER_H

#include "kmea.h"



static Qkmer* make_qkmer_from_str(const char* str, uint8_t length);
static Qkmer* qkmer_parse(const char* str);
static Qkmer* make_qkmer_from_kmer(Kmer* kmer);



#endif