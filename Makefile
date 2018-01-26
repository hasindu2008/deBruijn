CC     = gcc
LIBS = -Lhtslib -lz -lm -lbz2 -llzma -lpthread 
CFLAGS   = -g -Wall -O2
INC = -Ihtslib

all: debruijn.c common.c common.h
	$(CC) debruijn.c common.c htslib/libhts.a $(CLFLAGS) $(INC) $(LIBS) -o debruijn -g


clean: 
	rm debruijn
