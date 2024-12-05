#include <stdint.h>
#include <stdio.h>

void separate_kmer(uint8_t k, uint64_t kmer, uint64_t *ac, uint64_t *gt) {
    // Masks for extracting 2-bit pairs
    const uint64_t odd_position_pair_mask = 0x3333333333333333; // Binary: 00110011... (two 1s in each pair)
    const uint64_t even_position_pair_mask = 0xCCCCCCCCCCCCCCCC; // Binary: 11001100... (two 0s in each pair)
    const uint64_t zero_one_mask = 0x5555555555555555; // Binary: 01010101... 
    const uint64_t one_zero_mask = 0xAAAAAAAAAAAAAAAA; // Binary: 10101010... 

    print_binary(kmer);
    // Isolate T (11):
    uint64_t t_odd = kmer & odd_position_pair_mask;
    t_odd = t_odd & (t_odd >> 1); // t_odd now contains 01 where there are two 1s in the pair that is at an odd position
    uint64_t t_even = kmer & even_position_pair_mask;
    t_even = t_even & (t_even >> 1); // t_even now contains 10 where there are two 0s in the pair that is at an even position
    uint64_t t = t_odd | t_even; 
    printf("t: ");
    print_binary(t);

    // Isolate A (00): 
    uint64_t a_odd = ~kmer & odd_position_pair_mask; 
    uint64_t a_even = ~kmer & even_position_pair_mask;   
    uint64_t a = a_odd & (a_odd >> 1) | a_even & (a_even >> 1);
    a = a << 1; // to have 10 instead of 01
    printf("a: ");
    print_binary(a);

    // Isolate C (01):
    uint64_t c = kmer & zero_one_mask & ~t;
    printf("c: ");
    print_binary(c);

    // Isolate G (10):
    uint64_t g = (kmer & one_zero_mask) >> 1 & ~t;
    g = g << 1; // to have 10 instead of 01
    printf("g: ");
    print_binary(g);

    uint64_t length_mask = (1ULL << (k * 2)) - 1;

    *ac = (a | c) & length_mask; // Mask to get rid of the bits that are not part of the k-mer (is necessary here)
    *gt = (g | t) & length_mask; // Mask to get rid of the bits that are not part of the k-mer (should not shange anything here)

}
void print_binary(uint64_t value) {
    for (int i = 63; i >= 0; i--) {
        printf("%d", (value >> i) & 1);
        if (i % 8 == 0) printf(" "); // Add spaces for readability
        else if (i % 2 == 0) printf("."); // Add dots between 2-bit pairs
    }
    printf("\n");
}

int main() {
    uint64_t kmer = 0b111100011011; // Encoded "TTACGT"
    uint8_t k = 6;

	uint64_t complement = kmer ^ ((1 << k * 2) - 1);
    print_binary(complement);
    uint64_t reverse = 0;
    for (int i = 0; i < k; i++) {
        uint8_t byte = complement & 0b11;
        reverse <<= 2;
        reverse |= byte;
        complement >>= 2;
    }
    print_binary(reverse);
    return 0;
}
