CC=gcc
CFLAGS=-Wall
LDFLAGS=-lm
EXEC=volume
SRCS=$(wildcard *.c)
OBJS=$(patsubst %.c,%.o,$(SRCS))

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(LDFLAGS) $^ -o $@

%.o: %.c
	$(CC) $(CFLAGS) $< -c

clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC)
