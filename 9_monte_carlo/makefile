CC = gcc
CFLAGS = -Wall -lm -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

# Overall
all: A

clean:
	rm -f *.o *.out *.dat *.pdf mainA mainB mainC

# Object files
ode_rkf45.o: ../libs/ode_rkf45.c ../libs/ode_rkf45.h
	$(CC) $(CFLAGS) -c $< -o $@
monte_carlo.o: monte_carlo.c monte_carlo.h
adapt_2d.o: adapt_2d.c adapt_2d.h
# Task A
A: A.out A.pdf
A.out: mainA
	./mainA 1> A.out 2> A.dat
A.dat: A.out mainA
A.pdf: A.dat
	pyxplot A.ppl
mainA: mainA.c monte_carlo.o
# Task B
B: B.out
B.out: mainB
	./mainB > B.out
mainB: mainB.c monte_carlo.o adapt_2d.o

