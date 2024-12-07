DROP TABLE IF EXISTS kmers; 
DROP TABLE IF EXISTS DNAs;
DROP TABLE IF EXISTS qkmers;
DROP TABLE IF EXISTS large_table;
drop index if exists large_table_kmer_idx;
DROP EXTENSION IF EXISTS kmea;

CREATE EXTENSION kmea;
CREATE TABLE kmers(id INTEGER GENERATED ALWAYS AS IDENTITY PRIMARY KEY, kmer kmer); -- GENERATED ALWAYS AS IDENTITY = auto-increment for ID
CREATE TABLE DNAs(id INTEGER GENERATED ALWAYS AS IDENTITY PRIMARY KEY, dna DNA);
CREATE TABLE qkmers(id INTEGER GENERATED ALWAYS AS IDENTITY PRIMARY KEY, qkmer qkmer);

\copy DNAs (dna) FROM 'filtered_input.tsv' DELIMITER E'\t' CSV;

--SELECT * FROM DNAs;
SELECT COUNT(*) FROM DNAs;  --Ok with test "wc -l filtered_input.tsv"

SELECT id, length(dna) AS sequence_length FROM DNAs WHERE id <= 5;  --OK with check on fichier.fasta. We see we deleted sequence 4 beacause contained a N


SELECT COUNT(k.kmer) 
FROM generate_kmers((SELECT dna FROM DNAs WHERE id = 1), 5) AS k(kmer);
--Ok because Nb Kmer = L-k+1 | L = length of DNA sequence | k = length of kmer
-- Nb Kmer = 177-5+1 = 173


--INSERT INTO kmers (kmer)
INSERT INTO kmers (kmer)
SELECT k.kmer
FROM generate_kmers((SELECT dna FROM DNAs WHERE id = 1), 5) AS k(kmer);

SELECT * FROM kmers WHERE id <= 5;
SELECT * FROM kmers WHERE equals('ACGTT', kmer);
SELECT * FROM kmers WHERE kmer = 'ACGTA';   
SELECT * FROM kmers WHERE kmer ^@ 'ACGTA';  
SELECT * FROM kmers WHERE kmer ^@ 'ACGT';   --Ok because we have a kmer 'ACGTT' in kmers table

select COUNT(kmer) from kmers WHERE length(kmer)=5; --OK length of kmer = 5

select kmer as "Matches ACGNW" from kmers where 'ACGNW' @> kmer;    --OK because we see kmer 'ACGTT'


-- Create the table
CREATE TABLE large_table(
    id serial primary key, 
    kmer kmer
);

INSERT INTO large_table(kmer)
SELECT k.kmer
from generate_series(1, 200) as i, generate_kmers((SELECT dna FROM DNAs WHERE id = i), 5) as k(kmer);


create index large_table_kmer_idx on large_table using spgist(kmer spgist_kmer_ops);

SET enable_seqscan = on;
-- EXPLAIN ANALYSE SELECT count(*) FROM large_table where kmer ^@ 'AC';
SELECT count(*) FROM large_table where kmer ^@ 'AC';
SET enable_seqscan = off;

-- EXPLAIN ANALYSE SELECT count(*) FROM large_table where kmer ^@ 'AC';
SELECT count(*) FROM large_table where kmer ^@ 'AC';