CC = gcc
CFLAGS = -Wall -lm -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

# Overall
all: A B C

clean:
	rm -f *.o *.dat *.pdf mainA mainB mainC

# Object files
adapt_2d_speclim.o: adapt_2d_speclim.c adapt_2d_speclim.h
adapt_nd_speclim.o: adapt_nd_speclim.c adapt_nd_speclim.h
monte_carlo_speclim.o: monte_carlo_speclim.c monte_carlo_speclim.h

# Task A
A: A.out
A.out: mainA
	./mainA > $@
mainA: mainA.c adapt_2d_speclim.o

# Task B
B: B.out
B.out: mainB
	./mainB > $@
mainB: mainB.c monte_carlo_speclim.o adapt_nd_speclim.o

# Task C
C: C.out C_error.pdf C_calls.pdf
C.out: mainC
	./mainC > $@
mainC: mainC.c monte_carlo_speclim.o adapt_nd_speclim.o
C_error.pdf: C.out C.ppl
	pyxplot C.ppl
C_calls: C.out C.ppl
	pyxplot C.ppl

