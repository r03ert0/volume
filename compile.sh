if [ -f main.o ]; then rm main.o; fi
if [ -f Analyze.o ]; then rm Analyze.o; fi
if [ -f MGH.o ]; then rm MGH.o; fi
if [ -f Nifti.o ]; then rm Nifti.o; fi
if [ -f swap.o ]; then rm swap.o; fi

gcc -Wall -c main.c Analyze.c MGH.c Nifti.c swap.c
gcc -Wall -lm main.o Analyze.o MGH.o Nifti.o swap.o -o volume
