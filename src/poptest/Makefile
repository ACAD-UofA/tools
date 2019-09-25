TARGET=poptest.so
OBJS=poptest.o binom.o
HTSDIR=../htslib
BCFDIR=../bcftools
CFLAGS=-I$(BCFDIR) -I$(HTSDIR) -fPIC -Wall -O2 -ggdb3 -g -std=c99
LDFLAGS=-shared -L$(HTSDIR)
LIBS=-lm -lhts

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJS) -o $(TARGET) $(LIBS)

binom.o: binom.h
2poptest.o: binom.h

clean:
	rm -f $(TARGET) $(OBJS)
