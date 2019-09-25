TARGET=af_drift_test.so
OBJS=drift.o af_drift_test.o
HTSDIR=../htslib
BCFDIR=../bcftools
CFLAGS=-I$(BCFDIR) -I$(HTSDIR) -I./cephes -fPIC -Wall -O2 -ggdb3
LDFLAGS=-L./cephes -L$(HTSDIR)
LIBS=-lm -lcephes -lhts

all: $(TARGET) drift

$(TARGET): $(OBJS) cephes
	$(CC) $(CFLAGS) $(LDFLAGS) -shared $(OBJS) -o $(TARGET) $(LIBS)

drift: cephes
	$(CC) $(CFLAGS) -DTEST1 $(LDFLAGS) drift.c -o drift -lm -lcephes

.PHONY: cephes
cephes:
	$(MAKE) -C cephes

clean:
	rm -f $(TARGET) $(OBJS) drift
	$(MAKE) -C cephes clean
