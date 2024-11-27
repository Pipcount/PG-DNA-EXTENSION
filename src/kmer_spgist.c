#include "kmer.h"
#include "qkmer.h"
#include "access/spgist.h"

PG_MODULE_MAGIC;

typedef struct spgNodePtr
{
	Datum		d;
	int			i;
	int16		c;
} spgNodePtr;


/**
 * @brief Converts a K-mer to a string.
 * 
 * @param kmer The K-mer to convert.
 * @return The string representation of the K-mer.
 */
char* kmer_value_to_string(bytea* kmer) {
	//! elog(INFO, "kmer value (kmer_value_to_string): %lu", kmer -> value);
    uint8_t* data = (uint8_t*) VARDATA(kmer);
    uint8_t k = *data;
    elog(INFO, "kmer_value_to_string: k: %d", k);
    uint8_t first_byte_length = k % 4 == 0 ? 4 : k % 4;
	char* str = palloc0(k);
    data++;
    // Skip first (8 - ceil(k // 4)) bytes
    data += (8 - (int)ceil((double)k / 4));
    elog(INFO, "kmer_value_to_string: data: %x", *data);
	for (uint32_t i = 0; i < k; i += 4) {
        uint8_t current_byte = *data;
        data++;
        uint8_t number_of_nucleotides_to_read = i == 0 ? first_byte_length : 4; // read less nucleotides for the first byte

        for (int j = number_of_nucleotides_to_read - 1; j >= 0;  j--) {
            uint8_t nucleotide = (current_byte >> (j * 2)) & 0b11;                  // shift to the right by 6 bits for the first nucleotide, 4 bits for the second, etc.
            elog(INFO, "nucleotide: %d", nucleotide);
            str[i + j] = BINARY_TO_NUCLEOTIDE[nucleotide];
        }
    }
    str[k] = '\0';
    return str;
}

/**
 * @brief Finds the common prefix between two K-mers.
 * 
 * @param a The first K-mer value.
 * @param b The second K-mer value.
 * @param lena The k of the first K-mer.
 * @param lenb The k of the second K-mer.
 * @return The length of the common prefix.
 */
static int common_prefix(const uint8_t *a, const uint8_t *b, int lena, int lenb) {
	int i = 0;
    uint8_t a_val;
    uint8_t b_val;
    elog(INFO, "commonPrefix: lena: %d, lenb: %d", lena, lenb);
    elog(INFO, "commonPrefix: a: %d, b: %d", *a, *b);
	while (i < lena && i < lenb) {
		a++;
		b++;
		i++;
	}

	return i;
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

    cfg->prefixType = BYTEAOID;
    cfg->leafType = BYTEAOID;
    cfg->labelType = BYTEAOID;
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

    // Convert datum to k-mer (int64)
    Kmer *kmer = DatumGetKmerP(in->datum);
    elog(INFO, "kmer_spgist_choose: kmer->value: %s, kmer->k; %d", kmer_value_to_string(kmer), kmer->k);

    // Extract level-specific prefix (first N bits, where N depends on 'level')
    uint64_t current_prefix = (kmer->value >> (2 * (kmer->k - in->level - 1))) & 0b11; // Extract 2 bits for this level
    elog(INFO, "current_prefix: %d", current_prefix);
    // Handle case where current tuple has a prefix
    elog(INFO, "in->hasPrefix: %d", in->hasPrefix);
    if (in->hasPrefix) {
        uint64_t node_prefix = DatumGetUInt64(in->prefixDatum);

        if (node_prefix == current_prefix) {
            // Case 1: Prefix matches, descend into child
            out->resultType = spgMatchNode;
            out->result.matchNode.nodeN = 0;  // Assuming one node per prefix
            out->result.matchNode.levelAdd = 1; // Increment level by 1
            out->result.matchNode.restDatum = in->leafDatum; // Pass the same value down
            PG_RETURN_VOID();
        } else {
            // Case 3: Prefix doesn't match, split the node
            out->resultType = spgSplitTuple;
            out->result.splitTuple.prefixHasPrefix = true;
            out->result.splitTuple.prefixPrefixDatum = Int16GetDatum(current_prefix);
            out->result.splitTuple.prefixNNodes = in->nNodes + 1; // Add one new node
            out->result.splitTuple.prefixNodeLabels = (Datum *) palloc(sizeof(Datum) * out->result.splitTuple.prefixNNodes);

            // Copy existing node labels and add the new one
            for (int i = 0; i < in->nNodes; i++) {
                out->result.splitTuple.prefixNodeLabels[i] = in->nodeLabels[i];
            }
            out->result.splitTuple.prefixNodeLabels[in->nNodes] = Int16GetDatum(current_prefix); // Add new label

            out->result.splitTuple.childNodeN = in->nNodes; // This is the new node
            out->result.splitTuple.postfixHasPrefix = false; // No prefix for child node
            PG_RETURN_VOID();
        }
    }

    // Case 2: Add new node if no prefix or child matches
    out->resultType = spgAddNode;
    out->result.addNode.nodeLabel = Int64GetDatum(current_prefix);
    elog(INFO, "out->result.addNode.nodeLabel: %d", DatumGetInt64(out->result.addNode.nodeLabel));
    out->result.addNode.nodeN = in->nNodes; // Add at the end of the node list
    elog(INFO, "out->result.addNode.nodeN: %d", out->result.addNode.nodeN);
    PG_RETURN_VOID();
}

PG_FUNCTION_INFO_V1(kmer_spgist_picksplit);
Datum kmer_spgist_picksplit(PG_FUNCTION_ARGS) {
    elog(INFO, "kmer_spgist_picksplit");
    spgPickSplitIn *in = (spgPickSplitIn *) PG_GETARG_POINTER(0);
    spgPickSplitOut *out = (spgPickSplitOut *) PG_GETARG_POINTER(1);

    //// elog(INFO, "picksplit called: nTuples=%d, level=%d", in->nTuples, in->level);
    //// for (int i = 0; i < in->nTuples; i++) {
    ////     bytea *kmer = DatumGetByteaP(in->datums[i]);
    ////     uint8_t *data = (uint8_t *) VARDATA(kmer);
    ////     uint8_t k = *data;
    ////     elog(INFO, "k is: %d", k);
    ////     data++;
    ////     for (int j = 0; j < VARSIZE(kmer) - VARHDRSZ - 1; j++) {
    ////         elog(INFO, "j: %d, i: %d", j, *data);
    ////         data++;
    ////     }
    //// }

    bytea *kmer0 = DatumGetByteaP(in->datums[0]);
    uint8_t* data = (uint8_t*) VARDATA(kmer0);
    uint8_t k = *data;
    elog(INFO, "kmer0: %s", kmer_value_to_string(kmer0));
    int i, common_prefix_len;
    spgNodePtr* nodes;

    common_prefix_len = VARSIZE_ANY_EXHDR(kmer0) - 1;  // k value is stored in the first byte
    for (i = 1; i < in->nTuples && common_prefix_len > 0; i++) {
        bytea* kmeri = DatumGetByteaPP(in->datums[i]);
        int tempPrefixLen = common_prefix(VARDATA_ANY(kmer0),
									       VARDATA_ANY(kmeri),
									       VARSIZE_ANY_EXHDR(kmer0),
									       VARSIZE_ANY_EXHDR(kmeri));

        if (tempPrefixLen < common_prefix_len) {
            common_prefix_len = tempPrefixLen;
        }
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

/**
 * @brief Compress function to convert kmer to INT8OID.
 * 
 * @param fcinfo The function call information.
 * @return The compressed value.
 */
PG_FUNCTION_INFO_V1(kmer_spgist_compress);
Datum kmer_spgist_compress(PG_FUNCTION_ARGS) {
    Kmer *input = PG_GETARG_KMER_P(0);
    bytea *result = palloc0(VARHDRSZ + sizeof(uint64_t) + sizeof(uint8_t));        // + sizeof(uint8_t) for k value
    SET_VARSIZE(result, VARHDRSZ + sizeof(uint64_t) + sizeof(uint8_t));

    uint8_t* data_ptr = (uint8_t*) VARDATA(result);
    *data_ptr = input->k;
    data_ptr++;
    
    for (int8_t i = 56; i >= 0; i -= 8) {       // Start from 56 and decrement by 8 (8 bits per byte) (64 - 8 = 56)
        *data_ptr = (input->value >> i) & 0xFF;
        data_ptr++;
    }

    PG_RETURN_BYTEA_P(result);
}