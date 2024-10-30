#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "postgres.h"
#include "fmgr.h"

#include "libpq/pqformat.h"
#include "utils/fmgrprotos.h"

PG_MODULE_MAGIC;
typedef struct Kmer {
	uint64_t value;	// kmer value
	uint8_t k;  // kmer length
} Kmer;

#define DatumGetKmerP(X)  ((Kmer *) DatumGetPointer(X))
#define KmerPGetDatum(X)  PointerGetDatum(X)
#define PG_GETARG_KMER_P(n) DatumGetKmerP(PG_GETARG_DATUM(n))
#define PG_RETURN_KMER_P(x) return KmerPGetDatum(x)




// static inline kmer_t kmerFromStr(const char* str) {

// 	int str_len = strlen(str);
// 	if (str_len > 32) {
// 		perror("kmerFromStr: kmer is too long");
// 	}

// 	kmer_t kmer = 0;
	

// }
static Kmer* make_kmer(const char *str, uint8_t length) {
	Kmer* kmer = palloc0(sizeof(Kmer));
	kmer -> k = length;
	kmer -> value = 0;

	for (uint8_t i = 0; i < length; i++) {
		char c = str[i];
		if (c == 'A') {
			kmer -> value = (kmer -> value << 2) | 0b00;
		} else if (c == 'C') {
			kmer -> value = (kmer -> value << 2) | 0b01;
		} else if (c == 'G') {
			kmer -> value = (kmer -> value << 2) | 0b10;
		} else if (c == 'T') {
			kmer -> value = (kmer -> value << 2) | 0b11;
		} else {
			ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
	  	errmsg("invalid nucleotide")));
		}
	}
	return kmer;
}

static Kmer* kmer_parse(const char* str) {
	uint8_t length = strlen(str);
	if (length > 32) {
		ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
      	errmsg("kmer should not exceed 32 nucleotides")));
	}
	return make_kmer(str, length);
}

static char* kmer_value_to_string(Kmer* kmer) {

	char str[32];
	for (uint8_t i = 0; i < kmer -> k; i++) {
        uint8_t shift = (kmer -> k - i - 1) * 2;
        uint8_t nucleotide = (kmer -> value >> shift) & 0b11;

        if (nucleotide == 0b00) {
            str[i] = 'A';
        } else if (nucleotide == 0b01) {
            str[i] = 'C';
        } else if (nucleotide == 0b10) {
            str[i] = 'G';
        } else if (nucleotide == 0b11) {
            str[i] = 'T';
        }
    }
	str[kmer->k] = '\0';
	return psprintf("%s", str);
}
/* ************************************************************************** */

PG_FUNCTION_INFO_V1(kmer_in);
Datum kmer_in(PG_FUNCTION_ARGS) {
	char *str = PG_GETARG_CSTRING(0);
	PG_RETURN_KMER_P(kmer_parse(str));
}

PG_FUNCTION_INFO_V1(kmer_out);
Datum kmer_out(PG_FUNCTION_ARGS) {
	Kmer *kmer = PG_GETARG_KMER_P(0);
	char *str = kmer_value_to_string(kmer);
	PG_FREE_IF_COPY(kmer, 0);
	PG_RETURN_CSTRING(str);
}

PG_FUNCTION_INFO_V1(kmer_cast_from_text);
Datum
kmer_cast_from_text(PG_FUNCTION_ARGS)
{
  text *txt = PG_GETARG_TEXT_P(0);
  char *str = DatumGetCString(DirectFunctionCall1(textout,
               PointerGetDatum(txt)));
  PG_RETURN_KMER_P(kmer_parse(str));
}

PG_FUNCTION_INFO_V1(kmer_cast_to_text);
Datum
kmer_cast_to_text(PG_FUNCTION_ARGS)
{
  Kmer* kmer  = PG_GETARG_KMER_P(0);
  text* out = (text *)DirectFunctionCall1(textin,
            PointerGetDatum(kmer_value_to_string(kmer)));
  PG_FREE_IF_COPY(kmer, 0);
  PG_RETURN_TEXT_P(out);
}

