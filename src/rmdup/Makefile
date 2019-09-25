TARGET=rmdup
OBJS=rmdup.o
CFLAGS=-Wall -g -O2
LDFLAGS=
LIBS=

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $< -o $(TARGET) $(LDFLAGS) $(LIBS)

clean:
	rm -f $(OBJS) $(TARGET)
