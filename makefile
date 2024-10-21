EXTENSION   = kmea
MODULES     = kmea
DATA        = kmea--1.O.sql kmea.control

LDFLAGS=-lrt

PG_CONFIG ?= pg_config
PGXS = $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)
