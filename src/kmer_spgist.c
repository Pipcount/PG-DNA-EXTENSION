#include "kmer.h"
#include "qkmer.h"
#include "access/spgist.h"
#include "utils/pg_locale.h"
#include "utils/datum.h"

PG_MODULE_MAGIC;

#define EQUAL_STRATEGY_NUMBER 1
#define PREFIX_STRATEGY_NUMBER 2
#define QKMER_MATCHING_STRATEGY_NUMBER 3

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
    if (k > kmer->k) {
        ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE), errmsg("k cannot be greater than the K-mer size")));
    }
    Kmer* first_kmer = palloc0(sizeof(Kmer));
    first_kmer->k = k;
    first_kmer->value = kmer->value >> (2 * (kmer->k - k));
    return first_kmer;
}

/**
 * @brief Get the first k nucleotides of a QK-mer.
 * 
 * @param qkmer The QK-mer.
 * @param k The number of nucleotides to get.
 * @return The QK-mer with the first k nucleotides.
 */
Qkmer* get_first_k_nucleotides_qkmer(Qkmer* qkmer, uint8_t k) {
    if (k > qkmer->k) {
        ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE), errmsg("k cannot be greater than the QK-mer size")));
    }
    Qkmer* first_qkmer = palloc0(sizeof(Qkmer));
    first_qkmer->k = k;
    first_qkmer->ac = qkmer->ac >> (2 * (qkmer->k - k));
    first_qkmer->gt = qkmer->gt >> (2 * (qkmer->k - k));
    return first_qkmer;
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
 * @brief Compare the n first nucleotides of two K-mers.
 * 
 * @param kmer1 The first K-mer.
 * @param kmer2 The second K-mer.
 * @param n The number of nucleotides to compare.
 * @return The comparison result, -1 if kmer1 < kmer2, 0 if kmer1 == kmer2, 1 if kmer1 > kmer2.
 */
int compare_kmers(Kmer* kmer1, Kmer* kmer2, uint8_t n) {
    Kmer* first_kmer1 = get_first_k_nucleotides(kmer1, n);
    Kmer* first_kmer2 = get_first_k_nucleotides(kmer2, n);

    int result = 0;
    if (first_kmer1->value < first_kmer2->value) {
        result = -1;
    } else if (first_kmer1->value > first_kmer2->value) {
        result = 1;
    }
    pfree(first_kmer1);
    pfree(first_kmer2);
    return result;
}

/**
 * @brief Creates a Q-kmer from a K-mer.
 * 
 * @param kmer The K-mer to create the Q-kmer from.
 * @return A pointer to the created Q-kmer.
 */
static Qkmer* make_qkmer_from_kmer(Kmer* kmer) {
    Qkmer* qkmer = palloc0(sizeof(Qkmer));
    qkmer -> k = kmer -> k;

    // Masks for extracting 2-bit pairs
    const uint64_t odd_position_pair_mask = 0x3333333333333333; // Binary: 00110011... (two 1s in each pair)
    const uint64_t even_position_pair_mask = 0xCCCCCCCCCCCCCCCC; // Binary: 11001100... (two 0s in each pair)
    const uint64_t zero_one_mask = 0x5555555555555555; // Binary: 01010101... 
    const uint64_t one_zero_mask = 0xAAAAAAAAAAAAAAAA; // Binary: 10101010... 

    // Isolate T (11):
    uint64_t t_odd = kmer->value & odd_position_pair_mask;
    uint64_t t_even = kmer->value & even_position_pair_mask;
    uint64_t t = t_odd & (t_odd >> 1) | t_even & (t_even >> 1);

    // Isolate A (00): 
    uint64_t a_odd = ~kmer->value & odd_position_pair_mask; 
    uint64_t a_even = ~kmer->value & even_position_pair_mask;   
    uint64_t a = a_odd & (a_odd >> 1) | a_even & (a_even >> 1);
    a = a << 1;                                                     // to have 10 instead of 01

    // Isolate C (01):
    uint64_t c = kmer->value & zero_one_mask & ~t;

    // Isolate G (10):
    uint64_t g = (kmer->value & one_zero_mask) >> 1 & ~t;
    g = g << 1;                                                     // to have 10 instead of 01

    uint64_t length_mask = (1ULL << (kmer -> k * 2)) - 1;

    qkmer -> ac = (a | c) & length_mask; // Mask to get rid of the bits that are not part of the k-mer (is necessary here)
    qkmer -> gt = (g | t) & length_mask; // Mask to get rid of the bits that are not part of the k-mer (should not shange anything here)
    return qkmer;
}


/**
 * @brief Match a QK-mer with a K-mer.
 * 
 * @param qkmer The QK-mer.
 * @param kmer The K-mer.
 * @return true if the QK-mer matches the K-mer, false otherwise.
 */
bool qkmer_contains(Qkmer* qkmer, Kmer* kmer) {
    Qkmer* qkmer_from_kmer = make_qkmer_from_kmer(kmer);

    if (qkmer -> k != kmer -> k) {
        return false;
    }

    bool result = (qkmer_from_kmer -> ac & qkmer -> ac) == qkmer_from_kmer -> ac &&
                  (qkmer_from_kmer -> gt & qkmer -> gt) == qkmer_from_kmer -> gt;
    pfree(qkmer_from_kmer);
    return result;
}

/**
 * @brief Match a QK-mer with a K-mer up to a certain length.
 * 
 * @param qkmer The QK-mer.
 * @param kmer The K-mer.
 * @param n The number of nucleotides to compare.
 * @return true if the QK-mer matches the K-mer, false otherwise.
 */
bool qkmer_contains_n(Qkmer* qkmer, Kmer* kmer, uint8_t n) {
    Qkmer* first_qkmer = get_first_k_nucleotides_qkmer(qkmer, n);
    Kmer* first_kmer = get_first_k_nucleotides(kmer, n);

    bool result = qkmer_contains(first_qkmer, first_kmer);
    pfree(first_qkmer);
    pfree(first_kmer);
    return result;
}

/**
 * @brief Checks if a K-mer starts with a prefix.
 * 
 * @param kmer The K-mer to check.
 * @param prefix The prefix to check.
 * @return True if the K-mer starts with the prefix, false otherwise.
 */
static bool internal_kmer_startswith(Kmer* kmer, Kmer* prefix) {
	if (kmer -> k < prefix -> k) {
		return false;
	}
	uint64_t extracted_from_kmer = kmer -> value >> (kmer -> k - prefix -> k) * 2;
	return extracted_from_kmer == prefix -> value;
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
        if (first_non_common_nucleotide >= 0) {
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
    spgInnerConsistentIn *in = (spgInnerConsistentIn *) PG_GETARG_POINTER(0);
	spgInnerConsistentOut *out = (spgInnerConsistentOut *) PG_GETARG_POINTER(1);
    Kmer* prefix_kmer = NULL;
    Kmer* reconstructed_kmer;

    Kmer* reconstructed_value = DatumGetKmerP(in->reconstructedValue);
    Assert(reconstructed_value == NULL ? in->level == 0 : reconstructed_value->k == in->level);

    int max_reconstruction_length = in->level + 1;

    if (in->hasPrefix) {    // If we have a prefix
        prefix_kmer = DatumGetKmerP(in->prefixDatum);
        max_reconstruction_length += prefix_kmer->k;
    }
    reconstructed_kmer = palloc0(sizeof(Kmer));
    reconstructed_kmer->k = max_reconstruction_length - 1;
    reconstructed_kmer->value = 0;
    if (in->level) {        // If we are not at the root, we need to add the first nucleotides to the reconstructed K-mer
        Kmer* first_k_nucleotides = get_first_k_nucleotides(reconstructed_value, in->level);
        reconstructed_kmer->value = first_k_nucleotides->value;
        pfree(first_k_nucleotides);
    }
    /* 
     * Prefix_kmer could be null, if we have no prefix, so we first check that it's not
     * Then, if we have a prefix, we need to add it to the reconstructed K-mer
     */
    if (prefix_kmer && prefix_kmer->k) {
        reconstructed_kmer->value = (reconstructed_kmer->value << (2 * in->level)) | prefix_kmer->value;
    }

    /*
	 * Scan the child nodes.  For each one, complete the reconstructed value
	 * and see if it's consistent with the query.  If so, emit an entry into
	 * the output arrays.
     * https://github.com/postgres/postgres/blob/5d39becf8ba0080c98fee4b63575552f6800b012/src/backend/access/spgist/spgtextproc.c#L479
	 */
    out->nodeNumbers = (int *) palloc(sizeof(int) * in->nNodes);
	out->levelAdds = (int *) palloc(sizeof(int) * in->nNodes);
	out->reconstructedValues = (Datum *) palloc(sizeof(Datum) * in->nNodes);
	out->nNodes = 0;
    for (int i = 0; i < in->nNodes; i++) {
        int16 node_label = DatumGetInt16(in->nodeLabels[i]);
        Kmer* current_reconstructed_kmer_to_check = palloc0(sizeof(Kmer));
        *current_reconstructed_kmer_to_check = *reconstructed_kmer;
        bool result = true;
        if (node_label < 0) {
            current_reconstructed_kmer_to_check->k = max_reconstruction_length - 1;
        } else {
            current_reconstructed_kmer_to_check->k = max_reconstruction_length;
            current_reconstructed_kmer_to_check->value = (current_reconstructed_kmer_to_check->value << 2) | node_label;
        }
        for (int j = 0; j < in->nkeys; j++) {
            StrategyNumber strategy = in->scankeys[j].sk_strategy;
            
            if (strategy == QKMER_MATCHING_STRATEGY_NUMBER) {
                Qkmer* qkmer_in = DatumGetQkmerP(in->scankeys[j].sk_argument);
                // Compare with our function
                result = qkmer_contains_n(qkmer_in, current_reconstructed_kmer_to_check, Min(current_reconstructed_kmer_to_check->k, qkmer_in->k));
            } else {
                Kmer* kmer_in = DatumGetKmerP(in->scankeys[j].sk_argument);
                int compare_result = compare_kmers(kmer_in, current_reconstructed_kmer_to_check, Min(current_reconstructed_kmer_to_check->k, kmer_in->k));
                switch (strategy) {
                    case EQUAL_STRATEGY_NUMBER:
                        if (compare_result != 0 || kmer_in->k < current_reconstructed_kmer_to_check->k) {
                            result = false;
                        }
                        break;
                    case PREFIX_STRATEGY_NUMBER:
                        if (compare_result != 0) {
                            result = false;
                        }
                        break;
                    default:
                        ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE), errmsg("unrecognized strategy number: %d", strategy)));
                        break;
                }
                elog(INFO, "result: %d", result);
            }
            if (!result) {
                break;
            }
        }
        if (result) {
            out->nodeNumbers[out->nNodes] = i;
            out->levelAdds[out->nNodes] = current_reconstructed_kmer_to_check->k - in->level;
            current_reconstructed_kmer_to_check->k = max_reconstruction_length;
            out->reconstructedValues[out->nNodes] = datumCopy(KmerPGetDatum(current_reconstructed_kmer_to_check), false, sizeof(Kmer));
            out->nNodes++;
        }
    }
    PG_RETURN_VOID();
}

PG_FUNCTION_INFO_V1(kmer_spgist_leaf_consistent);
Datum
kmer_spgist_leaf_consistent(PG_FUNCTION_ARGS) {
    spgLeafConsistentIn *in = (spgLeafConsistentIn *) PG_GETARG_POINTER(0);
	spgLeafConsistentOut *out = (spgLeafConsistentOut *) PG_GETARG_POINTER(1);

    Kmer* full_kmer;
    out->recheck = false;

    Kmer* leaf_kmer = DatumGetKmerP(in->leafDatum);
    Kmer* reconstructed_value = DatumGetKmerP(in->reconstructedValue);

    int full_length = in->level + leaf_kmer->k;        // Full length of the K-mer
    if (leaf_kmer->k == 0 && in->level > 0) {
        full_kmer = reconstructed_value;
        out->leafValue = KmerPGetDatum(reconstructed_value);
    } else {
        full_kmer = palloc0(sizeof(Kmer));
        full_kmer->k = full_length;
        full_kmer->value = (reconstructed_value->value << (2 * leaf_kmer->k)) | leaf_kmer->value; // Combine the reconstructed value with the leaf value
        out->leafValue = KmerPGetDatum(full_kmer);
    }

    bool result = true;
    for (int j = 0; j < in->nkeys; j++) {
        StrategyNumber strategy = in->scankeys[j].sk_strategy;
        
        if (strategy == QKMER_MATCHING_STRATEGY_NUMBER) {
            Qkmer* qkmer_in = DatumGetQkmerP(in->scankeys[j].sk_argument);
            // We want to check if the QK-mer contains the full K-mer exactly (no truncation)
            result = qkmer_contains(qkmer_in, full_kmer);  
        } else if (strategy == EQUAL_STRATEGY_NUMBER) {
            Kmer* kmer_in = DatumGetKmerP(in->scankeys[j].sk_argument);
            // Check if the K-mers are equal
            result = KMER_EQUAL(kmer_in, full_kmer);
        } else if (strategy == PREFIX_STRATEGY_NUMBER) {
            Kmer* kmer_in = DatumGetKmerP(in->scankeys[j].sk_argument);
            // Check if the K-mer starts with the prefix
            result = internal_kmer_startswith(full_kmer, kmer_in);
        } else {
            ereport(ERROR, (errcode(ERRCODE_INVALID_PARAMETER_VALUE), errmsg("unrecognized strategy number: %d", strategy)));
            result = false;
        }

        if (!result) {
            break;
        }

    }
    
    PG_RETURN_BOOL(result);
}
