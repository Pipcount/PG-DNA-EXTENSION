DROP TABLE IF EXISTS kmers; 
DROP TABLE IF EXISTS DNAs;
DROP TABLE IF EXISTS qkmers;
DROP TABLE IF EXISTS large_table;
drop index if exists large_table_kmer_idx;

DROP EXTENSION IF EXISTS kmea;

CREATE EXTENSION kmea;
CREATE TABLE kmers(id integer, kmer kmer);
CREATE TABLE DNAs(id integer, dna DNA);
CREATE TABLE qkmers(id integer, qkmer qkmer);

INSERT INTO kmers VALUES 
(1, 'ACGTT'),
(2, 'ACTGTTTTTTT'),
(3, 'A'),
(4, 'ATCTACTCATCATACTGATCGATTAGat'),
(5, 'ACGAC'),
(6, 'ACGGG');

SELECT * FROM kmers;
select length(kmer), kmer from kmers;
SELECT * FROM kmers WHERE equals('ACGTT', kmer);
SELECT * FROM kmers WHERE kmer = 'ACGTA';
SELECT * FROM kmers WHERE kmer ^@ 'ACGTA';
SELECT * FROM kmers WHERE kmer ^@ 'ACGT';

INSERT INTO DNAs VALUES 
(1, 'ATCGG'),
(2, 'ATATATATATATATATA');

INSERT INTO DNAs VALUES
(3, dna('AAG'));

SELECT * FROM DNAs;

select length(dna), dna
from DNAs;

select k.kmer
from generate_kmers('ACGTACGT', 4) as  k(kmer)
WHERE k.kmer ^@ 'ACGT';

INSERT INTO qkmers VALUES 
(1, 'ACM'),
(2, 'GTNW'),
(3, 'ANMAAAAAAAAAA');

select * from qkmers;

select kmer as "Matches ACGNW"
from kmers
where 'ACGNW' @> kmer;

select length(qkmer), qkmer
from qkmers;


/* Test the kmer_count function */
select count(*)
from DNAs;

select k.kmer, count(*)
from generate_kmers('ACTGACTG', 4) as k(kmer)
group by k.kmer
order by count(*) desc;

with mykmers as (
    select k.kmer, count(*)
    from generate_kmers('ACTGACTG', 4) as k(kmer)
    group by k.kmer
)
select sum(count) as "Total kmers", count(*) as "Distinct kmers", count(*) filter (where count = 1) as "Unique kmers"
from mykmers;



-- Create the table
CREATE TABLE large_table(
    id serial primary key, 
    kmer kmer
);

-- Insert 1,000,000 rows into the table
INSERT INTO large_table (kmer)
SELECT k.kmer
FROM generate_kmers('CCACTAACAATGATCCACTAACAATGATCCACTAACAATGATCCACTAACAATGATCCACTAACAATGATCCACTAACAATGATCCACTAACAATGATCCACTAACAATGATCCACTAACAATGATCCACTAACAATGATCCACTAACAATGATCCACTAACAATGATCCACTAACAATGATCCACTAACAATGATCCACTAACAATGATCCACTAACAATGATCCACTAACAATGATCCACTAACAATGATCCACTAACAATGATCCACTAACAATGATCCACTAACAATGATCCACTAACAATGAT', 2) as k(kmer);

create index large_table_kmer_idx on large_table using spgist(kmer spgist_kmer_ops);

insert into large_table (kmer)
select k.kmer
from generate_kmers('CCACT', 2) as k(kmer);

EXPLAIN ANALYZE SELECT * FROM large_table WHERE kmer = 'AA';
-- EXPLAIN ANALYZE SELECT * FROM large_table WHERE kmer ^@ 'AC';