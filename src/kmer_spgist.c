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
 * @brief Get the first k nucleotides of a K-mer.
 * 
 * @param kmer The K-mer.
 * @param k The number of nucleotides to get.
 * @return The K-mer with the first k nucleotides.
 */
Kmer* get_first_k_nucleotides(Kmer* kmer, uint8_t k) {
    Kmer* first_kmer = palloc0(sizeof(Kmer));
    first_kmer->k = k;
    first_kmer->value = kmer->value >> (2 * (kmer->k - k));
    return first_kmer;
}

/**
 * @brief Get the last k nucleotides of a K-mer.
 * 
 * @param kmer The K-mer.
 * @param k The number of nucleotides to get.
 * @return The K-mer with the last k nucleotides.
 */
Kmer* get_last_k_nucleotides(Kmer* kmer, uint8_t k) {
    Kmer* last_kmer = palloc0(sizeof(Kmer));
    last_kmer->k = k;
    last_kmer->value = kmer->value & ((1 << (2 * k)) - 1);
    return last_kmer;
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
 * @brief Binary search an array of int16 datums for a match to c
 * On success, *i gets the match location; on failure, it gets where to insert
 * https://github.com/postgres/postgres/blob/master/src/backend/access/spgist/spgtextproc.c#L158
 * 
 * @param nodeLabels The array of int16 datums.
 * @param nNodes The number of nodes.
 * @param c The value to search for.
 * @param i The index of the value.
 * @return true if the value is found, false otherwise.
 */
static bool search_nucleotide(Datum *nodeLabels, int nNodes, int16 c, int *i) {
	int	StopLow = 0, StopHigh = nNodes;

	while (StopLow < StopHigh) {
		int	StopMiddle = (StopLow + StopHigh) >> 1;
		int16 middle = DatumGetInt16(nodeLabels[StopMiddle]);

		if (c < middle) {
			StopHigh = StopMiddle;
        } else if (c > middle) {
			StopLow = StopMiddle + 1;
        } else {
			*i = StopMiddle;
			return true;
		}
	}
	*i = StopHigh;
	return false;
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
    spgChooseIn *in = (spgChooseIn *) PG_GETARG_POINTER(0);
    spgChooseOut *out = (spgChooseOut *) PG_GETARG_POINTER(1);

    // Kmer that will be indexed
    Kmer *kmer_in = DatumGetKmerP(in->datum);

    uint8_t common_prefix_len = 0;
    int16 first_non_common_nucleotide = 0;
    
    if (in->hasPrefix) {
        // K-mer that is the prefix of the indexed K-mer
        Kmer* prefix_kmer = DatumGetKmerP(in->prefixDatum);

        // prefix_kmer at the current level
        Kmer* remaining_kmer_in = get_last_k_nucleotides(kmer_in, kmer_in->k - in->level);

        // common prefix length between the prefix K-mer and the indexed K-mer
        common_prefix_len = get_common_prefix_len(remaining_kmer_in, kmer_in);

        if (common_prefix_len == prefix_kmer->k) {
            if (common_prefix_len == kmer_in->k) {
                first_non_common_nucleotide = -1;
            } else {
                first_non_common_nucleotide = kmer_in->value >> (2 * (kmer_in->k - common_prefix_len - in->level - 1)) & 0b11;
            }
        } else {
            // split tuple
            out->resultType = spgSplitTuple;

            if (common_prefix_len == 0) {
                out->result.splitTuple.prefixHasPrefix = false;
            } else {
                out->result.splitTuple.prefixHasPrefix = true;
                Kmer* prefix_prefix_kmer = get_first_k_nucleotides(prefix_kmer, common_prefix_len);
                out->result.splitTuple.prefixPrefixDatum = KmerPGetDatum(prefix_prefix_kmer);
            }
            out->result.splitTuple.prefixNNodes = 1;
            out->result.splitTuple.prefixNodeLabels = (Datum *) palloc(sizeof(Datum));
            int16 nucleotide_after_prefix_prefix = prefix_kmer->value >> (2 * (prefix_kmer->k - common_prefix_len - 1)) & 0b11;
            out->result.splitTuple.prefixNodeLabels[0] = Int16GetDatum(nucleotide_after_prefix_prefix);
            out->result.splitTuple.childNodeN = 0;

            if (prefix_kmer->k - common_prefix_len == 1) {
                out->result.splitTuple.postfixHasPrefix = false;
            } else {
                out->result.splitTuple.postfixHasPrefix = true;
                Kmer* prefix_postfix_kmer = get_last_k_nucleotides(prefix_kmer, prefix_kmer->k - common_prefix_len - 1);
                out->result.splitTuple.postfixPrefixDatum = KmerPGetDatum(prefix_postfix_kmer);
            }
            PG_RETURN_VOID();
        }
    } else if (kmer_in->k > in->level) { 
        first_non_common_nucleotide = kmer_in->value >> (2 * (kmer_in->k - in->level - 1)) & 0b11;
    } else {
        first_non_common_nucleotide = -1;
    }
    int position = 0;
    if (search_nucleotide(in->nodeLabels, in->nNodes, first_non_common_nucleotide, &position)) {
       // Descent to existing node because it exists, at position "position"
        out->resultType = spgMatchNode;
        out->result.matchNode.nodeN = position;

        int level_add = common_prefix_len;
        if (first_non_common_nucleotide != -1) {
            level_add++;
        }
        out->result.matchNode.levelAdd = level_add;
        if (kmer_in->k - in->level - level_add > 0) {
            Kmer* remaining_kmer_in = get_last_k_nucleotides(kmer_in, kmer_in->k - in->level - level_add);
            out->result.matchNode.restDatum = KmerPGetDatum(remaining_kmer_in);
        } else {
            out->result.matchNode.restDatum = KmerPGetDatum(get_last_k_nucleotides(kmer_in, 0));
        }
    } else if (in->allTheSame) {
        // https://github.com/postgres/postgres/blob/master/src/backend/access/spgist/spgtextproc.c#L158
        out->resultType = spgSplitTuple;
		out->result.splitTuple.prefixHasPrefix = in->hasPrefix;
		out->result.splitTuple.prefixPrefixDatum = in->prefixDatum;
		out->result.splitTuple.prefixNNodes = 1;
		out->result.splitTuple.prefixNodeLabels = (Datum *) palloc(sizeof(Datum));
		out->result.splitTuple.prefixNodeLabels[0] = Int16GetDatum(-2);
		out->result.splitTuple.childNodeN = 0;
		out->result.splitTuple.postfixHasPrefix = false;
    } else {
		out->resultType = spgAddNode;
        out->result.addNode.nodeLabel = Int16GetDatum(first_non_common_nucleotide);
        out->result.addNode.nodeN = position;
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
    spgPickSplitIn *in = (spgPickSplitIn *) PG_GETARG_POINTER(0);
    spgPickSplitOut *out = (spgPickSplitOut *) PG_GETARG_POINTER(1);

    // for (int i = 0; i < in->nTuples; i++) {
    //     Kmer *kmer = DatumGetKmerP(in->datums[i]);
    //     elog(INFO, "kmer_spgist_picksplit: kmer->value: %s, kmer->k; %d", kmer_value_to_string(kmer), kmer->k);
    // }
    uint8_t common_prefix_len = get_common_prefix_len_array(in->datums, in->nTuples);

    if (common_prefix_len == 0) {
        out->hasPrefix = false;
    } else {
        out->hasPrefix = true;
        Kmer* kmer0 = DatumGetKmerP(in->datums[0]);
        Kmer* prefix_kmer = get_first_k_nucleotides(kmer0, common_prefix_len);
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
        if (common_prefix_len < kmeri->k) {
            Kmer* leaf_kmer = get_last_k_nucleotides(kmeri, kmeri->k - common_prefix_len - 1);
            leaf_datum = KmerPGetDatum(leaf_kmer);
        } else {
            Kmer* leaf_kmer = get_last_k_nucleotides(kmeri, 0);
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
