TARGET=vcfhet hom_windows
ZLIB=/opt/shared/zlib/1.2.8-gnu_4.8.0/lib
HTSDIR=../../../htslib
HTSFLAGS=-I$(HTSDIR)
HTSLIBS=-L$(HTSDIR) -L$(ZLIB) -lm -Wl,-static -lhts -lz -Wl,-Bdynamic -pthread
CFLAGS=-Wall -O2 -g $(HTSFLAGS)
LDLIBS=-lm

all: $(TARGET)

.SECONDEXPANSION:
vcfhet: $$@.o
	$(CC) $(CFLAGS) $^ -o $@ $(HTSLIBS) $(LDLIBS)

hom_windows: $$@.o
	$(CC) $(CFLAGS) $^ -o $@ $(HTSLIBS) $(LDLIBS)

clean:
	rm -f $(TARGET) $(addsuffix .o,$(TARGET))
