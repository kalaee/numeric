CC = gcc
CFLAGS = -Wall -lm -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

# Overall

all: A B

clean:
	rm -f *.o *.out mainA mainB_cyc mainB_row mainB_max

# Object files
jacobi.o: jacobi.c

# Task A
A: mainA
	./mainA > A.out

mainA: mainA.c jacobi.o jacobi.h

#Task B
B: B.out
mainB_row: mainB_row.c jacobi.o jacobi.h
mainB_max: mainB_max.c jacobi.o jacobi.h
mainB_cyc: mainB_cyc.c jacobi.o jacobi.h
B.out: mainB_row mainB_cyc mainB_max
	echo "Cyclic sweeps:" > $@
	for N in `seq 50 50 300`; do \
	echo "$$N"; \time -ao $@ -f "\t\ttime: %U" ./mainB_cyc $$N >> $@; \
	done
	echo "\nRow-wise:" >> B.out
	for N in `seq 50 50 300`; do \
	echo "$$N"; \time -ao $@ -f "\t\ttime %U" ./mainB_row $$N >> $@; \
	done
	echo "\nMaximum element of row:" >> B.out
	for N in `seq 50 50 300`; do \
	echo "$$N"; \time -ao $@ -f "\t\ttime %U" ./mainB_max $$N >> $@; \
	done

