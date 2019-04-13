CC=gcc
CFLAGS=-Wall
LDFLAGS=-lm
EXEC=volume

all: $(EXEC)

volume: main.o Analyze.o MGH.o Nifti.o
	$(CC) $(LDFLAGS) main.o Analyze.o MGH.o Nifti.o -o volume

main.o: main.c
	$(CC) $(CFLAGS) -c main.c

Analyze.o: Analyze.c
	$(CC) $(CFLAGS) -c Analyze.c

MGH.o: MGH.c
	$(CC) $(CFLAGS) -c MGH.c

Nifti.o: Nifti.c
	$(CC) $(CFLAGS) -c Nifti.c

clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC)
