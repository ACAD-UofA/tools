TARGET=dist
OBJS=dist.o
HTSDIR=/localscratch/Programs/htslib
ZLIB=/opt/shared/zlib/1.2.8-gnu_4.8.0
CFLAGS=-I$(HTSDIR) -I$(ZLIB)/include -fPIC -Wall -O2 -g #-O0 -ggdb3
LDFLAGS=-L$(HTSDIR) -L$(ZLIB)/lib

# statically link htslib (requiring zlib and pthreads)
LIBS=-Wl,-static -lhts -lz -Wl,-Bdynamic -pthread
# dynamically link
#LIBS=-lhts

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(LDFLAGS) $(OBJS) -o $(TARGET) $(LIBS)

clean:
	rm -f $(TARGET) $(OBJS)
