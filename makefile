EXTENSION   = kmea
MODULES     = src/kmer
DATA        = kmea--1.0.sql kmea.control

PG_CONFIG ?= pg_config
PGXS = $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)