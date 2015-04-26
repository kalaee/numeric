CC = gcc
CFLAGS = -Wall -lm -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

# All
all: A B C

clean:
	rm -f *.o *.out *.pdf mainA mainB mainC

# Object files
lspline.o: lspline.c
qspline.o: qspline.c
cspline.o: cspline.c

# Project A
A: A.dat plotA.ppl
	pyxplot plotA.ppl

A.dat: mainA
	./mainA

mainA: mainA.c qspline.h lspline.h lspline.o qspline.o

# Project B
B: B.dat plotB.ppl
	pyxplot plotB.ppl

B.dat: mainB
	./mainB

mainB: mainB.c qspline.h qspline.o

# Project B
C: C.dat plotC.ppl
	pyxplot plotC.ppl

C.dat: mainC
	./mainC

mainC: mainC.c cspline.h cspline.o
