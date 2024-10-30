CREATE OR REPLACE FUNCTION kmer_in(cstring)
RETURNS kmer
AS '$libdir/kmer'
LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION kmer_out(kmer)
RETURNS cstring
AS '$libdir/kmer'
LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION kmer_recv(internal)
RETURNS kmer
AS '$libdir/kmer'
LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION kmer_send(kmer)
RETURNS bytea
AS '$libdir/kmer'
LANGUAGE C IMMUTABLE STRICT;

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
LANGUAGE C IMMUTABLE STRICT;

CREATE OR REPLACE FUNCTION text(kmer)
RETURNS text
AS '$libdir/kmer', 'kmer_cast_to_text'
LANGUAGE C IMMUTABLE STRICT;

CREATE CAST (text as kmer) WITH FUNCTION kmer(text) AS IMPLICIT;
CREATE CAST (kmer as text) WITH FUNCTION text(kmer);
