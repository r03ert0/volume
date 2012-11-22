if [ -f main.o ]; then rm main.o; fi
if [ -f Analyze.o ]; then rm Analyze.o; fi
if [ -f MGH.o ]; then rm MGH.o; fi

gcc -Wall -c main.c Analyze.c MGH.c
gcc -Wall -lm main.o Analyze.o MGH.o -o volume
