-- -------------- --
-- Kmer data type --
-- -------------- --
CREATE OR REPLACE FUNCTION kmer_in(cstring)
RETURNS kmer
AS '$libdir/kmer'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION kmer_out(kmer)
RETURNS cstring
AS '$libdir/kmer'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION kmer_recv(internal)
RETURNS kmer
AS '$libdir/kmer'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION kmer_send(kmer)
RETURNS bytea
AS '$libdir/kmer'
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
AS '$libdir/kmer', 'kmer_cast_from_text'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION text(kmer)
RETURNS text
AS '$libdir/kmer', 'kmer_cast_to_text'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE CAST (text as kmer) WITH FUNCTION kmer(text) AS IMPLICIT;
CREATE CAST (kmer as text) WITH FUNCTION text(kmer);

CREATE OR REPLACE FUNCTION length(kmer)
RETURNS integer
AS '$libdir/kmer', 'kmer_length'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

-- OPERATORS
CREATE OR REPLACE FUNCTION equals(kmer, kmer)
RETURNS boolean
AS '$libdir/kmer', 'kmer_eq'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OPERATOR = (
	PROCEDURE = equals,
	LEFTARG = kmer,
	RIGHTARG = kmer,
	COMMUTATOR = =
);

CREATE OR REPLACE FUNCTION kmer_startswith(kmer, kmer)
RETURNS boolean
AS '$libdir/kmer', 'kmer_startswith'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OPERATOR ^@ (
	PROCEDURE = kmer_startswith,
	LEFTARG = kmer,
	RIGHTARG = kmer
);
	



-- -------------- --
-- DNA data type  --
-- -------------- --
CREATE OR REPLACE FUNCTION dna_in(cstring)
RETURNS DNA
AS '$libdir/dna'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION dna_out(DNA)
RETURNS cstring
AS '$libdir/dna'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION dna_recv(internal)
RETURNS DNA
AS '$libdir/dna'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION dna_send(DNA)
RETURNS bytea
AS '$libdir/dna'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE TYPE DNA (
    INPUT = dna_in,
    OUTPUT = dna_out,
    RECEIVE = dna_recv,
    SEND = dna_send
);

CREATE OR REPLACE FUNCTION DNA(text)
RETURNS DNA
AS '$libdir/dna', 'DNA_cast_from_text'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION text(DNA)
RETURNS text
AS '$libdir/dna', 'DNA_cast_to_text'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE CAST (text as DNA) WITH FUNCTION DNA(text) AS IMPLICIT;
CREATE CAST (DNA as text) WITH FUNCTION text(DNA);


CREATE OR REPLACE FUNCTION length(DNA)
RETURNS integer
AS '$libdir/dna', 'dna_length'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;


CREATE OR REPLACE FUNCTION generate_kmers(DNA, integer)
RETURNS SETOF kmer
AS '$libdir/dna', 'dna_generate_kmers'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

-- -------------- --
-- qkmer data type  --
-- -------------- --

CREATE OR REPLACE FUNCTION qkmer_in(cstring)
RETURNS qkmer
AS '$libdir/qkmer'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION qkmer_out(qkmer)
RETURNS cstring
AS '$libdir/qkmer'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION qkmer_recv(internal)
RETURNS qkmer
AS '$libdir/qkmer'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE OR REPLACE FUNCTION qkmer_send(qkmer)
RETURNS bytea
AS '$libdir/qkmer'
LANGUAGE C IMMUTABLE STRICT PARALLEL SAFE;

CREATE TYPE qkmer (
	INPUT = qkmer_in,
	OUTPUT = qkmer_out,
	RECEIVE = qkmer_recv,
	SEND = qkmer_send
);