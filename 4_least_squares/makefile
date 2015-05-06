CC = gcc
CFLAGS = -Wall -lm -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

# Overall

all: A B

clean:
	rm -f *.o *.out *.dat *.pdf mainA mainB

# Object files
qr.o: qr.c backsub.h
backsub.o: backsub.c
leasqr.o: leasqr.c backsub.h
jacobi.o: jacobi.c
svd.o: svd.c jacobi.h

# Task A
A: A.dat plotA.pdf

plotA.pdf: plotA.ppl A_points.dat A_fit.dat
	pyxplot plotA.ppl

A_fit.dat: mainA
	./mainA > /dev/null

A_points.dat: mainA
	./mainA > /dev/null

A.dat: mainA
	./mainA > A.out

mainA: mainA.c qr.o qr.h backsub.o leasqr.h leasqr.o

# Task B
B: B.dat plotB.pdf

plotB.pdf: plotB.ppl B_points.dat B_fit.dat
	pyxplot plotB.ppl

B_fit.dat: mainB
	./mainB > /dev/null

B_points.dat: mainB
	./mainB > /dev/null

B.dat: mainB
	./mainB > B.out

mainB: mainB.c qr.o qr.h backsub.o svd.o svd.h jacobi.h jacobi.o
