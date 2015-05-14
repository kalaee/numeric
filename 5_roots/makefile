CC = gcc
CFLAGS = -Wall -lm -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

# Overall
all: A B C

clean:
	rm -f *.o *.out mainA mainB mainC

# Object files
backsub.o: ../2_lineq/backsub.c ../2_lineq/backsub.h
	$(CC) $(CLFAGS) -c $< -o $@
givens.o: ../2_lineq/givens.c ../2_lineq/givens.h
	$(CC) $(CFLAGS) -c $< -o $@
vector.o: ../libs/vector.c ../libs/vector.h
	$(CC) $(CFLAGS) -c $< -o $@
newton.o: newton.c newton.h

# Task A
A: mainA
	./mainA 1> A.out
mainA: mainA.c givens.o backsub.o newton.o vector.o

# Task B
B: mainB
	./mainB 1> B.out
mainB: mainB.c givens.o backsub.o newton.o vector.o

# Task C
C: mainC
	./mainC 1> C.out
mainC: mainC.c givens.o backsub.o newton.o vector.o

