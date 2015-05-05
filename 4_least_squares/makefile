CC = gcc
CFLAGS = -Wall -lm -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

# Overall

all: A

clean:
	rm -f *.o *.out mainA

# Object files
qr.o: qr.c backsub.h
backsub.o: backsub.c

# Task A
A: mainA
	./mainA > A.out

mainA: mainA.c qr.o qr.h backsub.o
