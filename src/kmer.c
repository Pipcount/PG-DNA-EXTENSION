#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "postgres.h"
#include "fmgr.h"

#include "libpq/pqformat.h"
#include "utils/fmgrprotos.h"
#include "utils/elog.h"

PG_MODULE_MAGIC;
typedef struct Kmer {
	uint64_t value;	// kmer value
	uint8_t k;  // kmer length
} Kmer;

#define DatumGetKmerP(X)  ((Kmer *) DatumGetPointer(X))
#define KmerPGetDatum(X)  PointerGetDatum(X)
#define PG_GETARG_KMER_P(n) DatumGetKmerP(PG_GETARG_DATUM(n))
#define PG_RETURN_KMER_P(x) return KmerPGetDatum(x)

static Kmer* make_kmer(const char *str, uint8_t length) {
	Kmer* kmer = palloc0(sizeof(Kmer));
	kmer -> k = length;
	kmer -> value = 0;

	for (uint8_t i = 0; i < length; i++) {
		char c = str[i];
		switch (c) {
			case 'A': kmer -> value = (kmer -> value << 2) | 0b00; break;
			case 'C': kmer -> value = (kmer -> value << 2) | 0b01; break;
			case 'G': kmer -> value = (kmer -> value << 2) | 0b10; break;
			case 'T': kmer -> value = (kmer -> value << 2) | 0b11; break;
			default:
				ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
				errmsg("invalid nucleotide")));
				break;
		}
	}
	//! elog(INFO, "kmer value after make_kmer: %lu", kmer -> value);
	return kmer;
}

static Kmer* kmer_parse(const char* str) {
	uint8_t length = strlen(str);
	//! elog(INFO, "kmer length (kmer_parse): %d", length);
	if (length > 32) {
		ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
      	errmsg("kmer should not exceed 32 nucleotides")));
	} else if (length == 0) {
		ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
	  	errmsg("kmer should not be empty")));
	}
	return make_kmer(str, length);
}

static char* kmer_value_to_string(Kmer* kmer) {
	char str[32];
	//! elog(INFO, "kmer value (kmer_value_to_string): %lu", kmer -> value);
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
	Kmer* kmer = kmer_parse(str);
	//! elog(INFO, "kmer value (kmer_in): %lu", kmer -> value);
	PG_RETURN_KMER_P(kmer);
}

PG_FUNCTION_INFO_V1(kmer_out);
Datum kmer_out(PG_FUNCTION_ARGS) {
	Kmer *kmer = PG_GETARG_KMER_P(0);
	//! elog(INFO, "kmer value (kmer_out): %lu, size %d", kmer -> value, kmer -> k);
	char *str = kmer_value_to_string(kmer);
	//! elog(INFO, "kmer string (kmer_out): %sblablabla", str);
	PG_FREE_IF_COPY(kmer, 0);
	PG_RETURN_CSTRING(str);
}

PG_FUNCTION_INFO_V1(kmer_recv);
Datum kmer_recv(PG_FUNCTION_ARGS) {
	//! elog(INFO, "kmer_recv");
	StringInfo buf = (StringInfo) PG_GETARG_POINTER(0);
	Kmer *kmer = palloc(sizeof(Kmer));
	kmer -> value = pq_getmsgint64(buf);
	kmer -> k = pq_getmsgint(buf, 8);
	PG_RETURN_KMER_P(kmer);
}

PG_FUNCTION_INFO_V1(kmer_send);
Datum kmer_send(PG_FUNCTION_ARGS) {
	//! elog(INFO, "kmer_send");
	Kmer *kmer = PG_GETARG_KMER_P(0);
	StringInfoData buf;
	pq_begintypsend(&buf);
	pq_sendint64(&buf, kmer -> value);
	pq_sendint8(&buf, kmer -> k);
	PG_FREE_IF_COPY(kmer, 0);
	PG_RETURN_BYTEA_P(pq_endtypsend(&buf));
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

PG_FUNCTION_INFO_V1(kmer_length);
Datum
kmer_length(PG_FUNCTION_ARGS)
{
  Kmer* kmer  = PG_GETARG_KMER_P(0);
  PG_RETURN_CHAR(kmer -> k);
}