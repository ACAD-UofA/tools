#TARGET=svg
TARGET=parsecols
CFLAGS=-Wall -g -O0
LDLIBS=

$(TARGET): $(TARGET).o svgfig.o
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

clean:
	rm -f $(TARGET) *.o
