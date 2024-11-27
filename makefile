EXTENSION   = kmea
MODULES     = src/kmer src/dna src/qkmer src/kmer_spgist
DATA        = kmea--1.0.sql kmea.control

PG_CONFIG ?= pg_config
PGXS = $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)
