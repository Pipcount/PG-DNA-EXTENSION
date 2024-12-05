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

COPY DNAs (dna)
FROM '/mnt/c/Users/theod/OneDrive/Documents/ULB/Ma1/INFOH417-DatabaseSystemsArchitecture/Projet/PG-DNA-EXTENSION/filtered_input.tsv'
DELIMITER E'\t' --Indique que le délimiteur est une tabulation (\t), qui est couramment utilisée dans les fichiers TSV (Tab-Separated Values).
CSV;    --or TSV

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