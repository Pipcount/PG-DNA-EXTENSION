-- -------------- --
-- Kmer data type --
-- -------------- --
CREATE OR REPLACE FUNCTION kmer_in(cstring)
RETURNS kmer
AS '$libdir/kmea'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION kmer_out(kmer)
RETURNS cstring
AS '$libdir/kmea'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION kmer_recv(internal)
RETURNS kmer
AS '$libdir/kmea'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION kmer_send(kmer)
RETURNS bytea
AS '$libdir/kmea'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE TYPE kmer (
	INPUT = kmer_in,
	OUTPUT = kmer_out,
	RECEIVE = kmer_recv,
	SEND = kmer_send,
	INTERNALLENGTH = 9
);

CREATE OR REPLACE FUNCTION kmer(text)
RETURNS kmer
AS '$libdir/kmea', 'kmer_cast_from_text'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION text(kmer)
RETURNS text
AS '$libdir/kmea', 'kmer_cast_to_text'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE CAST (text as kmer) WITH FUNCTION kmer(text) AS IMPLICIT;
CREATE CAST (kmer as text) WITH FUNCTION text(kmer);

CREATE OR REPLACE FUNCTION length(kmer)
RETURNS integer
AS '$libdir/kmea', 'kmer_length'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION canonical(kmer)
RETURNS kmer
AS '$libdir/kmea', 'kmer_canonical'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

-- OPERATORS
CREATE OR REPLACE FUNCTION equals(kmer, kmer)
RETURNS boolean
AS '$libdir/kmea', 'kmer_eq'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OPERATOR = (
	PROCEDURE = equals,
	LEFTARG = kmer,
	RIGHTARG = kmer,
	COMMUTATOR = =
);

CREATE OR REPLACE FUNCTION startswith(prefix kmer, kmer kmer)
RETURNS boolean
AS '$libdir/kmea', 'kmer_startswith'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION startswith_inv(kmer kmer, prefix kmer)
RETURNS boolean
AS '$libdir/kmea', 'kmer_startswith_inv'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OPERATOR ^@ (
	PROCEDURE = startswith_inv,
	LEFTARG = kmer,
	RIGHTARG = kmer
);
	



-- -------------- --
-- DNA data type  --
-- -------------- --
CREATE OR REPLACE FUNCTION dna_in(cstring)
RETURNS DNA
AS '$libdir/kmea'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION dna_out(DNA)
RETURNS cstring
AS '$libdir/kmea'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION dna_recv(internal)
RETURNS DNA
AS '$libdir/kmea'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION dna_send(DNA)
RETURNS bytea
AS '$libdir/kmea'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE TYPE DNA (
    INPUT = dna_in,
    OUTPUT = dna_out,
    RECEIVE = dna_recv,
    SEND = dna_send
);

CREATE OR REPLACE FUNCTION DNA(text)
RETURNS DNA
AS '$libdir/kmea', 'DNA_cast_from_text'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION text(DNA)
RETURNS text
AS '$libdir/kmea', 'DNA_cast_to_text'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE CAST (text as DNA) WITH FUNCTION DNA(text) AS IMPLICIT;
CREATE CAST (DNA as text) WITH FUNCTION text(DNA);


CREATE OR REPLACE FUNCTION length(DNA)
RETURNS integer
AS '$libdir/kmea', 'dna_length'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;


CREATE OR REPLACE FUNCTION generate_kmers(DNA, integer)
RETURNS SETOF kmer
AS '$libdir/kmea', 'dna_generate_kmers'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

-- -------------- --
-- qkmer data type  --
-- -------------- --

CREATE OR REPLACE FUNCTION qkmer_in(cstring)
RETURNS qkmer
AS '$libdir/kmea'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION qkmer_out(qkmer)
RETURNS cstring
AS '$libdir/kmea'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION qkmer_recv(internal)
RETURNS qkmer
AS '$libdir/kmea'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION qkmer_send(qkmer)
RETURNS bytea
AS '$libdir/kmea'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE TYPE qkmer (
	INPUT = qkmer_in,
	OUTPUT = qkmer_out,
	RECEIVE = qkmer_recv,
	SEND = qkmer_send,
	INTERNALLENGTH = 17
);

CREATE OR REPLACE FUNCTION qkmer(text)
RETURNS qkmer
AS '$libdir/kmea', 'qkmer_cast_from_text'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION text(qkmer)
RETURNS text
AS '$libdir/kmea', 'qkmer_cast_to_text'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE CAST (text as qkmer) WITH FUNCTION qkmer(text) AS IMPLICIT;
CREATE CAST (qkmer as text) WITH FUNCTION text(qkmer);

CREATE OR REPLACE FUNCTION contains(qkmer, kmer)
RETURNS boolean
AS '$libdir/kmea', 'qkmer_contains'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OPERATOR @> (
	PROCEDURE = contains,
	LEFTARG = qkmer,
	RIGHTARG = kmer
);

CREATE OR REPLACE FUNCTION length(qkmer)
RETURNS integer
AS '$libdir/kmea', 'qkmer_length'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;


-- ------------------ --
-- Kmer Hash opclass  --
-- ------------------ --

CREATE OR REPLACE FUNCTION kmer_hash(kmer)
RETURNS integer
AS '$libdir/kmea', 'kmer_hash'
LANGUAGE C IMMUTABLE;


CREATE OPERATOR CLASS hash_kmer_ops
DEFAULT FOR TYPE kmer USING hash
AS
        OPERATOR        1       =  ,
		FUNCTION 	  	1       kmer_hash(kmer);


-- ------------------- --
-- Kmer SP-GiST index  --
-- ------------------- --


CREATE OR REPLACE FUNCTION kmer_spgist_config(internal, internal)
RETURNS void
AS '$libdir/kmea', 'kmer_spgist_config'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION kmer_spgist_choose(internal, internal)
RETURNS void
AS '$libdir/kmea', 'kmer_spgist_choose'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION kmer_spgist_picksplit(internal, internal)
RETURNS void
AS '$libdir/kmea', 'kmer_spgist_picksplit'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION kmer_spgist_inner_consistent(internal, internal)
RETURNS void
AS '$libdir/kmea', 'kmer_spgist_inner_consistent'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION kmer_spgist_leaf_consistent(internal, internal)
RETURNS boolean
AS '$libdir/kmea', 'kmer_spgist_leaf_consistent'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;



CREATE OPERATOR CLASS spgist_kmer_ops
DEFAULT FOR TYPE kmer USING spgist
AS
    OPERATOR    1   = (kmer, kmer) ,
    OPERATOR    2   ^@(kmer, kmer) ,
    OPERATOR    3   @>(qkmer, kmer) ,
    FUNCTION    1   kmer_spgist_config(internal, internal),
    FUNCTION    2   kmer_spgist_choose(internal, internal),
    FUNCTION    3   kmer_spgist_picksplit(internal, internal),
    FUNCTION    4   kmer_spgist_inner_consistent(internal, internal),
    FUNCTION    5   kmer_spgist_leaf_consistent(internal, internal);