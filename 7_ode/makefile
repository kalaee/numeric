CC = gcc
CFLAGS = -Wall -lm -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

# Overall
all: A

clean:
	rm -f *.o *.out *.pdf mainA mainB mainC

# Object files
ode_rkf45.o: ode_rkf45.c ode_rkf45.h
ode_driver.o: ode_driver.c ode_driver.c

# Task A
A: A.pdf

A.pdf: A.out A.ppl
	pyxplot A.ppl

A.out: mainA
	./mainA > A.out
mainA: mainA.c ode_driver.o ode_rkf45.o

