#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))

typedef struct Kmer {
	uint64_t value;	// kmer value
	uint8_t k;  // kmer length
} Kmer;

typedef int64_t DNA;


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

// void print_binary(Kmer* kmer) {
//     for (int i = 8*sizeof(kmer->value); i > 0 ; i--){
//         int k = kmer->value >> (i - 1);
//         if (k & 1) {
//             printf("1");
//         } else {
//             printf("0");
//         }
//     }
//     printf("\n");
// }

void print_binary(DNA* Dna) {
    for (int i = 8*sizeof(Dna); i > 0 ; i--){
        int k = *Dna >> (i - 1);
        if (k & 1) {
            printf("1");
        } else {
            printf("0");
        }
    }
    printf("\n");
}


static void store_last_byte_mask(DNA* dna, uint32_t length) {
    uint8_t* data_ptr = malloc(sizeof(uint8_t));
    uint8_t last_byte_mask = 0b11111111;
    for (uint32_t i = 0; i < length % 4; i++) {
        last_byte_mask = (last_byte_mask << 2) | 0b11;
    }

    *data_ptr = last_byte_mask;
    printf("%d\n", *data_ptr);
    free(data_ptr);
}

static DNA* make_dna(const char* str, uint32_t length) {
    DNA* dna = malloc(sizeof(dna));
    for (uint32_t i = 0; i < length; i += 4) {
        uint8_t current_byte = 0b00000000;
        for (char j = 0; j < 4; j++) {
            char c = str[i + j];
            if (i + j >= length) {
                current_byte = (current_byte << 2) | 0b00;
                continue;
            }
            if (c == 'A') {
                current_byte = (current_byte << 2) | 0b00;
            } else if (c == 'C') {
                current_byte = (current_byte << 2) | 0b01;
            } else if (c == 'G') {
                current_byte = (current_byte << 2) | 0b10;
            } else if (c == 'T') {
                current_byte = (current_byte << 2) | 0b11;
            } else {
                printf("invalid nucleotide\n");
            }
        }
        *dna = *dna << 8 | current_byte;
    }
    return dna; 
}

uint8_t create_last_byte_mask(uint8_t len_remainder) {
    uint8_t mask = (1 << (len_remainder * 2)) - 1;
    return mask << (8 - (len_remainder * 2));
}

void dna_generate_kmers() {
    char* str = 'ATGCATGC';
    DNA* dna = make_dna(str, strlen(str));
    uint8_t kmer_length = 6;
    uint32_t length = strlen(dna);
    
    uint8_t* byte_ptr = *dna;
    uint8_t* start_ptr = *byte_ptr;
    uint8_t nucleotide_ctr = 0;
    
    while (byte_ptr - start_ptr < (length - kmer_length) / 4) {
        Kmer* kmer = malloc(sizeof(Kmer));
        kmer->value = 0;
        kmer->k = kmer_length;
        


    }

}

bool check_pairs(uint64_t num, int k) {

    // XOR pairs of bits - this trick will set the bits to 01 if the pair is 10 or 01
    uint64_t xor_pairs = (num ^ (num >> 1)) & 0x5555555555555555ULL;
    printf("After XOR: ");
    print_binary(&xor_pairs);

    // Create a mask with 01010101...01 (k times)
    uint64_t pairs_of_01 = 0x5555555555555555ULL >> (64 - 2*k);
    printf("Pairs of 01 of lenght k: ");
    print_binary(&pairs_of_01);

    return xor_pairs == pairs_of_01;
}

int main() {
    // uint8_t length = 4;
    // char* kmer_str = generate_kmer(length);
    // printf("kmer: %s\n", kmer_str);
    // Kmer* kmer = make_kmer(kmer_str, length);

    // print_binary(kmer);
    // print_kmer(kmer);
    // free(kmer);
    // free(kmer_str);

    // char* dna_seq = "TCGATC";
    // DNA* dna = make_dna(dna_seq, strlen(dna_seq));
    // print_binary(dna);
    // free(dna);

    //printf("Mask for 6 bits: %02X\n", create_last_byte_mask(3));  // 11111100 -> FC
    //printf("Mask for 2 bits: %02X\n", create_last_byte_mask(1));  // 11000000 -> C0
    //printf("Mask for 8 bits: %02X\n", create_last_byte_mask(4));  // 11111111 -> FF
    //printf("Mask for 8 bits: %02X\n", create_last_byte_mask(2));  // 11111111 -> FF
    uint64_t num = 0ULL;
    num |= 0b0001100101010101011010101010100101010101101010101010101010101010;
    printf("initial number: ");
    print_binary(&num);
    // true = 1, false = 0
    bool result = check_pairs(num, 32);
    printf("Result: %s\n", result ? "true" : "false");

    return 0;
}