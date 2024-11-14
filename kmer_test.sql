DROP TABLE IF EXISTS kmers; 
DROP TABLE IF EXISTS DNAs;
DROP EXTENSION IF EXISTS kmea;

CREATE EXTENSION kmea;
CREATE TABLE kmers(id integer, kmer kmer);
CREATE TABLE DNAs(id integer, dna DNA);

INSERT INTO kmers VALUES 
(1, 'ACGTT'),
(2, 'ACTGTTTTTTT'),
(3, 'A'),
(4, 'ATCTACTCATCATACTGATCGATTAGat');

SELECT * FROM kmers;
select length(kmer), kmer from kmers;
SELECT * FROM kmers WHERE equals('ACGTT', kmer);
SELECT * FROM kmers WHERE kmer = 'ACGTA';
SELECT * FROM kmers WHERE kmer ^@ 'ACGTA';
SELECT * FROM kmers WHERE kmer ^@ 'ACGT';

INSERT INTO DNAs VALUES 
(1, 'ATCGG'),
(2, 'ATATATATATATATATA');

SELECT * FROM DNAs;

select length(dna), dna
from DNAs;

select k.kmer
from generate_kmers('ACGTACGT', 4) as  k(kmer)
WHERE k.kmer ^@ 'ACGT';

