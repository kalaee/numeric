CC = gcc
CFLAGS = -Wall -lm -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

# Overall
all: A B C

clean:
	rm -f *.o *.dat *.out *.pdf mainA mainB mainC

# Object files
vector.o: ../libs/vector.c ../libs/vector.h
	$(CC) $(CLFAGS) -c $< -o $@
ode_rkf45.o: ode_rkf45.c ode_rkf45.h
ode_rk5.o: ode_rk5.c ode_rk5.h
ode_driver.o: ode_driver.c ode_driver.c
ode_integrate.o: ode_integrate.c ode_integrate.h

# Task A
A: A.pdf

A.pdf: A.out A.ppl
	pyxplot A.ppl

A.out: mainA
	./mainA > A.dat
mainA: mainA.c ode_driver.o ode_rkf45.o

# Task B
B: B.pdf B.out

B.pdf: B.out B.ppl
	pyxplot B.ppl
B.out: mainB
	./mainB 2> B.out 1> B.dat

mainB: mainB.c ode_driver.o ode_rkf45.o ode_rk5.o vector.o

# Task C
C: C.out

C.out: mainC
	./mainC > C.out

mainC: mainC.c ode_driver.o ode_rkf45.o ode_integrate.o
