EXTENSION   = kmea
MODULE_big  = $(EXTENSION)

objdir = bin
srcdir = src

OBJS_C  = kmer.o dna.o qkmer.o kmer_spgist.o
OBJS   = $(addprefix src/, $(OBJS_C))

INCS   = kmer.h dna.h qkmer.h kmea.h

DATA        = kmea--1.0.sql kmea.control

PG_CONFIG = pg_config
PGXS = $(shell $(PG_CONFIG) --pgxs)
include $(PGXS)

$(OBJS): $(addprefix src/, $(INCS))

ifdef VPATH
all: vpath-mkdirs
.PHONY: vpath-mkdirs
$(OBJS): | vpath-mkdirs

vpath-mkdirs:
	$(MKDIR_P) $(objdir)
endif