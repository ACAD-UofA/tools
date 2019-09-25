TARGET?=mn
OBJS=$(TARGET).o kmath.o
HTSDIR=../htslib
CFLAGS=-I$(HTSDIR) -O2 -Wall -g
LDLIBS=-L$(HTSDIR) -lm -lhts -pthread

$(TARGET): $(OBJS)

clean:
	rm -f $(TARGET) $(OBJS)
