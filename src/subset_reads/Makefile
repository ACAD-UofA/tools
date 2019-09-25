TARGET=subset_reads
CFLAGS=-Wall -O2 -g
LDLIBS=-lm

$(TARGET): $(TARGET).o kmath.o
	$(CC) $(CFLAGS) $^ -o $(TARGET) $(LDLIBS)

clean:
	rm $(TARGET) *.o
