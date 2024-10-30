#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

typedef struct Kmer {
	uint64_t value;	// kmer value
	uint8_t k;  // kmer length
} Kmer;


static Kmer* make_kmer(const char *str, uint8_t length) {
	Kmer* kmer = malloc(sizeof(Kmer));
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
            printf("invalid nucleotide");
		}
	}
	return kmer;
}

void print_kmer(Kmer* kmer) {
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
    printf("%s\n", str);
}

static char* generate_kmer(uint8_t length) {
    char* kmer = malloc(length + 1);
    for (uint8_t i = 0; i < length; i++) {
        uint8_t nucleotide = rand() % 4;
        if (nucleotide == 0) {
            kmer[i] = 'A';
        } else if (nucleotide == 1) {
            kmer[i] = 'C';
        } else if (nucleotide == 2) {
            kmer[i] = 'G';
        } else if (nucleotide == 3) {
            kmer[i] = 'T';
        }
    }
    kmer[length] = '\0';
    return kmer;
}

void print_binary(Kmer* kmer) {
    for (int i = 8*sizeof(kmer->value); i > 0 ; i--){
        int k = kmer->value >> (i - 1);
        if (k & 1) {
            printf("1");
        } else {
            printf("0");
        }
    }
    printf("\n");
}

int main() {
    uint8_t length = 32;
    char* kmer_str = generate_kmer(length);
    printf("kmer: %s\n", kmer_str);
    Kmer* kmer = make_kmer(kmer_str, length);

    print_binary(kmer);
    print_kmer(kmer);
    free(kmer);
    free(kmer_str);
    return 0;
}