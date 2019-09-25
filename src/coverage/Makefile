TARGET=coverage
OBJS=$(TARGET).o
HTSDIR=../../../htslib
ZLIB=/opt/shared/zlib/1.2.8-gnu_4.8.0/lib
CFLAGS=-I$(HTSDIR) -Wall -O2 -g
LDLIBS=-L$(HTSDIR) -L$(ZLIB) -lm -Wl,-static -lhts -lz -Wl,-Bdynamic -pthread

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

clean:
	rm -f $(TARGET) *.o
