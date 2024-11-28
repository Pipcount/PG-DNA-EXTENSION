#include "kmer.h"
#include "qkmer.h"
#include "access/spgist.h"

PG_MODULE_MAGIC;

typedef struct KmerNodePtr
{
	Datum		kmer_datum;
	int			position;
	int16		first_non_common_nucleotide;
} KmerNodePtr;


/**
 * @brief Converts a K-mer to a string.
 * 
 * @param kmer The K-mer to convert.
 * @return The string representation of the K-mer.
 */
char* kmer_value_to_string(Kmer* kmer) {
	char str[32];
	//! elog(INFO, "kmer value (kmer_value_to_string): %lu", kmer -> value);
	for (uint8_t i = 0; i < kmer -> k; i++) {
        uint8_t shift = (kmer -> k - i - 1) * 2;
        uint8_t nucleotide = (kmer -> value >> shift) & 0b11;

		str[i] = BINARY_TO_NUCLEOTIDE[nucleotide];
    }
	str[kmer->k] = '\0';
	return psprintf("%s", str);
}

/**
 * @brief Function to get the common prefix length of 2 K-mers.
 * 
 * @param kmer1 The first K-mer.
 * @param kmer2 The second K-mer.
 * @return The common prefix length.
 */
uint8_t get_common_prefix_len(Kmer* kmer1, Kmer* kmer2) {
    uint64_t kmer1_value = kmer1->value; // Copy the value to avoid modifying the original
    uint64_t kmer2_value = kmer2->value; // Copy the value to avoid modifying the original

    int size_diff = kmer1->k - kmer2->k;
    uint8_t prefix_len;

    /*
        * If the K-mers have different sizes, we need to align them by shifting the larger one.
        * We then compare the values to find the common prefix length.
    */
    if (size_diff > 0) {
        prefix_len = kmer2->k;
        kmer1_value >>= 2 * size_diff;
    } else if (size_diff < 0) {
        prefix_len = kmer1->k;
        kmer2_value >>= 2 * -size_diff;
    } else {
        prefix_len = kmer1->k;
    }
    uint64_t xor = kmer1_value ^ kmer2_value;
    while (xor != 0) {                            // if prefix_len is 0, xor will be 0
        xor >>= 2;
        prefix_len--;
    }
    return prefix_len;
}


/**
 * @brief Function to get the common prefix length of a set of K-mers.
 * 
 * @param datums The K-mers to compare.
 * @param nTuples The number of K-mers.
 * @return The common prefix length.
 */
uint8_t get_common_prefix_len_array(Datum *datums, int nTuples) {
    Kmer* first_kmer = DatumGetKmerP(datums[0]);
    uint8_t common_prefix_len = first_kmer->k;
    for (int i = 1; i < nTuples && common_prefix_len > 0; i++) {
        Kmer* kmer = DatumGetKmerP(datums[i]);

        uint8_t prefix_len = get_common_prefix_len(first_kmer, kmer);

        if (prefix_len < common_prefix_len) {
            common_prefix_len = prefix_len;
        }
    }
    return common_prefix_len;
}

        

/**
 * @brief Postgres function that defines the SP-GiST configuration for K-mers.
 * 
 * @param cfg The SP-GiST configuration.
 * @return void
 */
PG_FUNCTION_INFO_V1(kmer_spgist_config);
Datum kmer_spgist_config(PG_FUNCTION_ARGS) {
    spgConfigOut *cfg = (spgConfigOut *) PG_GETARG_POINTER(1);
    spgConfigIn *in = (spgConfigIn *) PG_GETARG_POINTER(0);
    elog(INFO, "kmer_spgist_config");

    cfg->prefixType = in->attType;
    cfg->labelType = INT2OID; // 2 bits for the nucleotide
    // leaftype not initialized see https://www.postgresql.org/docs/current/spgist.html#SPGIST-BUILTIN-OPCLASSES-TABLE
    cfg->canReturnData = true;
    cfg->longValuesOK = false;

    PG_RETURN_VOID();
}

/**
 * @brief Postgres function to choose the best K-mer to split a page.
 * 
 * @param fcinfo The function call information.
 * @return void
 */
PG_FUNCTION_INFO_V1(kmer_spgist_choose);
Datum
kmer_spgist_choose(PG_FUNCTION_ARGS)
{
    elog(INFO, "kmer_spgist_choose");
    spgChooseIn *in = (spgChooseIn *) PG_GETARG_POINTER(0);
    spgChooseOut *out = (spgChooseOut *) PG_GETARG_POINTER(1);
    Kmer *kmer_in = DatumGetKmerP(in->datum);
    uint64_t prefix_kmer = 0;
    uint8_t prefix_size = 0;
    uint8_t common_prefix_len = 0;

    
    if(in ->hasPrefix){
        Kmer *kmer = DatumGetKmerP(in->prefixDatum);
        prefix_kmer = kmer->value;
        prefix_size = kmer->k;

        common_prefix_len = get_common_prefix_len(kmer_in, kmer);
    }
    PG_RETURN_VOID();
}

static int compare_nodes(const void *a, const void *b) {
    const KmerNodePtr *node_a = (const KmerNodePtr *) a;
    const KmerNodePtr *node_b = (const KmerNodePtr *) b;
    return (int32) node_a->first_non_common_nucleotide - (int32) node_b->first_non_common_nucleotide;
}
    

PG_FUNCTION_INFO_V1(kmer_spgist_picksplit);
Datum kmer_spgist_picksplit(PG_FUNCTION_ARGS) {
    elog(INFO, "kmer_spgist_picksplit");
    spgPickSplitIn *in = (spgPickSplitIn *) PG_GETARG_POINTER(0);
    spgPickSplitOut *out = (spgPickSplitOut *) PG_GETARG_POINTER(1);

    // for (int i = 0; i < in->nTuples; i++) {
    //     Kmer *kmer = DatumGetKmerP(in->datums[i]);
    //     elog(INFO, "kmer_spgist_picksplit: kmer->value: %s, kmer->k; %d", kmer_value_to_string(kmer), kmer->k);
    // }
    uint8_t common_prefix_len = get_common_prefix_len_array(in->datums, in->nTuples);
    elog(INFO, "common_prefix_len: %d", common_prefix_len);

    if (common_prefix_len == 0) {
        out->hasPrefix = false;
    } else {
        out->hasPrefix = true;
        Kmer* kmer0 = DatumGetKmerP(in->datums[0]);
        Kmer* prefix_kmer = palloc0(sizeof(Kmer));
        prefix_kmer->k = common_prefix_len;
        prefix_kmer->value = kmer0->value >> (2 * (kmer0->k - common_prefix_len));
        out->prefixDatum = KmerPGetDatum(prefix_kmer);
    }

    // Extract the first non-common nucleotide for each K-mer (Node Label)
    KmerNodePtr* nodes = (KmerNodePtr *) palloc(sizeof(KmerNodePtr) * in->nTuples);

    for (int i = 0; i < in->nTuples; i++) {
        Kmer* kmer = DatumGetKmerP(in->datums[i]);
        if (common_prefix_len == kmer->k) {
            nodes[i].first_non_common_nucleotide = -1;
        } else {
            nodes[i].first_non_common_nucleotide = kmer->value >> (2 * (kmer->k - common_prefix_len - 1)) & 0b11;
        }
        nodes[i].kmer_datum = in->datums[i];
        nodes[i].position = i;
    }

    qsort(nodes, in->nTuples, sizeof(*nodes), compare_nodes);

    out->nNodes =  0;
    out->nodeLabels = (Datum *) palloc(sizeof(Datum) * in->nTuples);
    out->mapTuplesToNodes = (int *) palloc(sizeof(int) * in->nTuples);
    out->leafTupleDatums = (Datum *) palloc(sizeof(Datum) * in->nTuples);

    for (int i = 0; i < in->nTuples; i++) {
        Kmer* kmeri = DatumGetKmerP(nodes[i].kmer_datum);
        Datum leaf_datum;

        if (i == 0 || nodes[i].first_non_common_nucleotide != nodes[i - 1].first_non_common_nucleotide) {
            out->nodeLabels[out->nNodes] = Int16GetDatum(nodes[i].first_non_common_nucleotide);
            out->nNodes++;
        }
        Kmer* leaf_kmer = palloc0(sizeof(Kmer));
        if (common_prefix_len < kmeri->k) {
            leaf_kmer->k = kmeri->k - common_prefix_len - 1;
            leaf_kmer->value = kmeri->value & ((1 << (2 * leaf_kmer->k)) - 1); 
            leaf_datum = KmerPGetDatum(leaf_kmer);
        } else {
            leaf_kmer->k = 0;
            leaf_kmer->value = 0;
            leaf_datum = KmerPGetDatum(leaf_kmer);
        }

        out->leafTupleDatums[nodes[i].position] = leaf_datum;
        out->mapTuplesToNodes[nodes[i].position] = out->nNodes - 1;
    }

    PG_RETURN_VOID();
}


PG_FUNCTION_INFO_V1(kmer_spgist_inner_consistent);
Datum
kmer_spgist_inner_consistent(PG_FUNCTION_ARGS) {
    PG_RETURN_VOID();
}

PG_FUNCTION_INFO_V1(kmer_spgist_leaf_consistent);
Datum
kmer_spgist_leaf_consistent(PG_FUNCTION_ARGS) {
    PG_RETURN_BOOL(1);
}
