\set db_name 'kmea'

DROP DATABASE IF EXISTS :db_name;
CREATE DATABASE :db_name;
\c :db_name;

CREATE EXTENSION kmea;
CREATE TABLE kmers(id serial primary key, kmer kmer); 
CREATE TABLE DNAS(id serial primary key, dna DNA);

\COPY DNAS (dna) FROM 'filtered_input.tsv';

\set nb_sequences 50000
-- Fill the kmers table with kmers generated FROM the DNA sequences
WITH small_table AS (
    SELECT * 
    FROM DNAS 
    WHERE id <= :nb_sequences      -- 100k rows, gives a reasonable time for testing
)
INSERT INTO kmers(kmer)
SELECT k.kmer
FROM small_table, LATERAL generate_kmers(dna, 30) AS k(kmer);


-- Add a few small kmers because it will be easier to test the equality operator
WITH small_table AS (
    SELECT * 
    FROM DNAS 
    WHERE id <= (:nb_sequences/80)
)
INSERT INTO kmers(kmer)
SELECT k.kmer
FROM small_table, LATERAL generate_kmers(dna, 5) AS k(kmer);


-- Create an index on the kmers table
CREATE INDEX kmer_idx ON kmers USING spgist(kmer spgist_kmer_ops);


-- Test the generate_kmers function
-- For reference, print the DNA sequence FROM which the kmers are generated (here we know that the first sequence is big enough to generate 10 kmers)
SELECT dna AS "DNA sequence"
FROM DNAS
WHERE id = 1;

SELECT id, kmer AS "Generated kmer"
FROM kmers
WHERE id <= 10;


-- Test the equality operator and function
SELECT count(*) AS "Amount that equals ACTGC"
FROM kmers WHERE equals('ACTGC', kmer);

SELECT count(*) AS "Amount that equals ACTGC" 
FROM kmers WHERE kmer = 'ACTGC';


-- Test the startswith operator and function
SELECT count(*) AS "Amount that starts with ACTG" 
FROM kmers WHERE startswith('ACTG', kmer);

SELECT count(*) AS "Amount that starts with ACTG" 
FROM kmers WHERE kmer ^@ 'ACTG';


-- Test the contains operator and function
SELECT count(*) AS "Amount that matches HAWKNNSAYBBDHVDDNNNNDRVWMDDHDS" 
FROM kmers WHERE contains('HAWKNNSAYBBDHVDDNNNNDRVWMDDHDS', kmer);

SELECT count(*) AS "Amount that matches HAWKNNSAYBBDHVDDNNNNDRVWMDDHDS" 
FROM kmers WHERE 'HAWKNNSAYBBDHVDDNNNNDRVWMDDHDS' @> kmer;


-- Test the counting support
SELECT kmer, count(*) AS "Amount of occurrences of the k-mers"
FROM kmers
GROUP BY kmer
ORDER BY count(*) DESC
LIMIT 10;

WITH counted AS(
    SELECT kmer, count(*)
    FROM kmers
    GROUP BY kmer
)
SELECT sum(count) AS "Total count",
count(*) AS "Distinct count",
count(*) FILTER (WHERE count = 1) AS "Unique count"
FROM counted;


-- Test the index scan vs seq scan
SET enable_seqscan = on;
EXPLAIN ANALYZE SELECT count(*) FROM kmers WHERE kmer ^@ 'ACTGCA';
SELECT count(*) as "Amount that starts with ACTGCA using SEQ SCAN"
FROM kmers
WHERE kmer ^@ 'ACTGCA';
SET enable_seqscan = off;
EXPLAIN ANALYZE SELECT count(*) FROM kmers WHERE kmer ^@ 'ACTGCA';
SELECT count(*) as "Amount that starts with ACTGCA using INDEX SCAN"
FROM kmers
WHERE kmer ^@ 'ACTGCA';