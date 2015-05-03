CC = gcc
CFLAGS = -Wall -lm -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

# Overall

all: A B C time.out

clean:
	rm -f *.o *.out mainA mainB mainC

# Object files
qr.o: qr.c backsub.h
givens.o: givens.c backsub.h
cholesky.o: cholesky.c backsub.h
backsub.o: backsub.c

# Task A
A: mainA
	./mainA > A.out

mainA: mainA.c qr.o qr.h backsub.o

# Task B
B: mainB
	./mainB > B.out

mainB: mainB.c givens.o givens.h backsub.o

time.out: mainA mainB
	echo 'Comparison of runtime for algorithms in A and B' > time.out
	echo 'Both solve exactly the same problems and' >> time.out
	echo 'reserve the same amount of memory.' >> time.out
	echo 'Note that they both call the same GSL LU decomp' >> time.out
	echo 'to evaluate their results.' >> time.out
	echo 'QR with Gram-Schmidt orthogonalisation:' >> time.out
	\time -f "%U sec" -ao time.out ./mainA > /dev/null
	echo 'QR with Givens rotations' >> time.out
	\time -f "%U sec" -ao time.out ./mainB > /dev/null

# Task C
C: mainC
	./mainC > C.out

mainC: mainC.c cholesky.o cholesky.h backsub.o
