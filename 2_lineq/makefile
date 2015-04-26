CC = gcc
CFLAGS = -Wall -lm -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

# Overall

all: A B

clean:
	rm -f *.o *.out mainA mainB mainC

time: mainA mainB
	{ \time ./mainA ; } 2> times.out

# General
qr.o: qr.c backsub.h
givens.o: givens.c backsub.h
cholesky.o: cholesky.c backsub.h

# Task A
A: mainA
	./mainA > A.out
	cat A.out

mainA: mainA.c qr.o qr.h

# Task B
B: mainB
	./mainB > B.out
	cat B.out

mainB: mainB.c givens.o givens.h

# Task C
C: mainC
	./mainC > C.out
	cat C.out

mainC: mainC.c cholesky.o cholesky.h
