CC = gcc
CFLAGS = -Wall -lm -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

# Overall

all: A

clean:
	rm -f *.o *.out mainA

# Object files
jacobi.o: jacobi.c

# Task A
A: mainA
	./mainA > A.out

mainA: mainA.c jacobi.o jacobi.h
