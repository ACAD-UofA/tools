HTSDIR=../htslib
CFLAGS=-I$(HTSDIR) -O2 -Wall -g
LDLIBS=-L$(HTSDIR) -lm -lhts -pthread

all: count_snp_patterns count_snp_patterns2 hybrid hpower

hybrid: kmath.o
	$(CC) $(CFLAGS) -D_HYBRID_MAIN -lm -pthread kmath.o hybrid.c -o $@

hpower: kmath.o hybrid.o
	$(CC) $(CFLAGS) -D_POWER_MAIN -lm -pthread kmath.o hybrid.c -o $@

count_snp_patterns: count_snp_patterns.o
	$(CC) $(CFLAGS) -L$(HTSDIR) $^ -o $@ -Wl,-static -lhts -lz -Wl,-Bdynamic -pthread
count_snp_patterns2: count_snp_patterns2.o
	$(CC) $(CFLAGS) -L$(HTSDIR) $^ -o $@ -Wl,-static -lhts -lz -Wl,-Bdynamic -pthread

clean:
	rm -f hybrid hpower count_snp_patterns count_snp_patterns2 *.o
