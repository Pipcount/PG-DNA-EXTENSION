#include "qkmer.h"

PG_MODULE_MAGIC;


// A C G T W S M K R Y B D H V N
static Qkmer* make_qkmer_from_str(const char* str, uint8_t length) {
    Qkmer* qkmer = palloc0(sizeof(Qkmer));
    qkmer -> k = length;
    qkmer -> ac = 0;
    qkmer -> gt = 0;

    for (uint8_t i = 0; i < length; i++) {
        char c = str[i];
        switch (toupper(c)) {
            case 'A': qkmer -> ac = (qkmer -> ac << 2) | 0b10; qkmer -> gt = (qkmer -> gt << 2) | 0b00; break;
            case 'C': qkmer -> ac = (qkmer -> ac << 2) | 0b01; qkmer -> gt = (qkmer -> gt << 2) | 0b00; break;
            case 'G': qkmer -> ac = (qkmer -> ac << 2) | 0b00; qkmer -> gt = (qkmer -> gt << 2) | 0b10; break;
            case 'T': qkmer -> ac = (qkmer -> ac << 2) | 0b00; qkmer -> gt = (qkmer -> gt << 2) | 0b01; break;
            case 'W': qkmer -> ac = (qkmer -> ac << 2) | 0b10; qkmer -> gt = (qkmer -> gt << 2) | 0b01; break;
            case 'S': qkmer -> ac = (qkmer -> ac << 2) | 0b01; qkmer -> gt = (qkmer -> gt << 2) | 0b10; break;
            case 'M': qkmer -> ac = (qkmer -> ac << 2) | 0b11; qkmer -> gt = (qkmer -> gt << 2) | 0b00; break;
            case 'K': qkmer -> ac = (qkmer -> ac << 2) | 0b00; qkmer -> gt = (qkmer -> gt << 2) | 0b11; break;
            case 'R': qkmer -> ac = (qkmer -> ac << 2) | 0b10; qkmer -> gt = (qkmer -> gt << 2) | 0b10; break;
            case 'Y': qkmer -> ac = (qkmer -> ac << 2) | 0b01; qkmer -> gt = (qkmer -> gt << 2) | 0b01; break;
            case 'B': qkmer -> ac = (qkmer -> ac << 2) | 0b01; qkmer -> gt = (qkmer -> gt << 2) | 0b11; break;
            case 'D': qkmer -> ac = (qkmer -> ac << 2) | 0b10; qkmer -> gt = (qkmer -> gt << 2) | 0b11; break;
            case 'H': qkmer -> ac = (qkmer -> ac << 2) | 0b11; qkmer -> gt = (qkmer -> gt << 2) | 0b01; break;
            case 'V': qkmer -> ac = (qkmer -> ac << 2) | 0b11; qkmer -> gt = (qkmer -> gt << 2) | 0b10; break;
            case 'N': qkmer -> ac = (qkmer -> ac << 2) | 0b11; qkmer -> gt = (qkmer -> gt << 2) | 0b11; break;
            default:
                ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
                errmsg("invalid nucleotide")));
                break;
        }
    }
    return qkmer;
}


static Qkmer* qkmer_parse(const char* str) {
    uint8_t length = strlen(str);

	if (length > 32) {
		ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
      	errmsg("qkmer should not exceed 32 nucleotides")));
	} else if (length == 0) {
		ereport(ERROR, (errcode(ERRCODE_INVALID_TEXT_REPRESENTATION),
	  	errmsg("qkmer should not be empty")));
	}
	return make_qkmer_from_str(str, length);
}
/* ************************************************************************** */

PG_FUNCTION_INFO_V1(qkmer_in);
Datum qkmer_in(PG_FUNCTION_ARGS) {
    char* str = PG_GETARG_CSTRING(0);
    Qkmer* qkmer = qkmer_parse(str);
    PG_RETURN_QKMER_P(qkmer);
}

PG_FUNCTION_INFO_V1(qkmer_out);
Datum qkmer_out(PG_FUNCTION_ARGS) {

}

PG_FUNCTION_INFO_V1(qkmer_send);
Datum qkmer_send(PG_FUNCTION_ARGS) {

}

PG_FUNCTION_INFO_V1(qkmer_recv);
Datum qkmer_recv(PG_FUNCTION_ARGS) {

}