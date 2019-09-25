TARGET?=methcounts
CFLAGS=-Wall -O2 -g
LDLIBS=-lm

$(TARGET): $(TARGET).o
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

clean:
	rm -f $(TARGET) $(addsuffix .o,$(TARGET))
