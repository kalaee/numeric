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
roots_newton.o: ../libs/roots_newton.c ../libs/roots_newton.h
	$(CC) $(CFLAGS) -c $< -o $@
min_newton.o: min_newton.c min_newton.h
min_newton_sr1.o: min_newton_sr1.c min_newton_sr1.h
simplex.o: simplex.c simplex.h

# Task A
A: mainA
	./mainA > A.out
mainA: mainA.c givens.o backsub.o min_newton.o vector.o

# Task B
B: mainB
	./mainB > B.out
mainB: mainB.c backsub.o givens.o min_newton.o min_newton_sr1.o vector.o roots_newton.o

#Task C
C: mainC
	./mainC > C.out
mainC: mainC.c simplex.o
