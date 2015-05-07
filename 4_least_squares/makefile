CC = gcc
CFLAGS = -Wall -lm -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

# Overall
all: A B C

clean:
	rm -f *.o *.out *.dat *.pdf mainA mainB mainC

# Object files
qr.o: ../2_lineq/qr.c ../2_lineq/qr.h
	$(CC) $(CLFAGS) -c $< -o $@
backsub.o: ../2_lineq/backsub.c ../2_lineq/backsub.h
	$(CC) $(CLFAGS) -c $< -o $@
jacobi.o: ../3_eigen/jacobi.c ../3_eigen/jacobi.h
	$(CC) $(CLFAGS) -c $< -o $@
svd.o: svd.c svd.h
leasqr.o: leasqr.c leasqr.h

# Task A
A: A.out plotA.pdf
plotA.pdf: plotA.ppl A.dat
	pyxplot plotA.ppl
A.dat: mainA
	./mainA 2> A.dat
A.out: mainA
	./mainA 1> A.out 2> A.dat
mainA: mainA.c qr.o backsub.o leasqr.o

# Task B
B: B.out plotB.pdf
plotB.pdf: plotB.ppl B.dat
	pyxplot plotB.ppl
B.dat: mainB
	./mainB 2> B.dat
B.out: mainB
	./mainB 1> B.out 2> B.dat
mainB: mainB.c svd.o jacobi.o

# Task C
C: C.out plotC.pdf
plotC.pdf: plotC.ppl C.dat
	pyxplot plotC.ppl
C.dat: mainC
	./mainC 2> C.dat
C.out: mainC
	./mainC 1> C.out 2> C.dat
mainC: mainC.c qr.o backsub.o leasqr.o
