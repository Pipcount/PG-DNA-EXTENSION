DROP TABLE IF EXISTS t;
DROP EXTENSION IF EXISTS kmea;

CREATE EXTENSION kmea;
CREATE TABLE t(id integer, kmer kmer);

INSERT INTO t VALUES 
(1, 'ACGTT'),
(2, 'ACTGTTTTTTT'),
(3, 'A'),
(4, 'ATCTACTCATCATACTGATCGATTAGCTGATG');

SELECT * FROM t;

select length(kmer), kmer from t;