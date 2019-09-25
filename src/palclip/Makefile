TARGET=palclip
OBJS=palclip.o
CFLAGS=-Wall -g -O2
LDFLAGS=
LIBS=-Wl,-Bstatic -lz -Wl,-Bdynamic

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $< -o $(TARGET) $(LDFLAGS) $(LIBS)

clean:
	rm -f $(OBJS) $(TARGET)
