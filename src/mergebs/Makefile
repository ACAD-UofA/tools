PROG=mergehbs
OBJS=mergehbs.o fit_lognorm.o kmath.o fold.o
CFLAGS=-Wall -g -Os
LDFLAGS=$(CFLAGS)
LDLIBS=-lm -lz

all: mergehbs simhbs

mergehbs: mergehbs.o fit_lognorm.o kmath.o fold.o

simhbs: simhbs.o fold.o kmath.o

test_fit: test_fit.o fold.o kmath.o fit_lognorm.o

clean:
	rm -f $(PROG) $(OBJS)
