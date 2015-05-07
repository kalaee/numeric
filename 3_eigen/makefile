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
mainB_row: mainB_row.c jacobi.o jacobi.h
mainB_max: mainB_max.c jacobi.o jacobi.h
mainB_cyc: mainB_cyc.c jacobi.o jacobi.h

B: mainB_row mainB_cyc mainB_max
	echo "Cyclic sweeps:" > B.out
	\time -f "\t\ttime: %U seconds" -ao B.out ./mainB_cyc 10 >> B.out
	\time -f "\t\ttime: %U seconds" -ao B.out ./mainB_cyc 50 >> B.out
	\time -f "\t\ttime: %U seconds" -ao B.out ./mainB_cyc 100 >> B.out
	\time -f "\t\ttime: %U seconds" -ao B.out ./mainB_cyc 150 >> B.out
	echo "\nRow-wise:" >> B.out
	\time -f "\t\ttime: %U seconds" -ao B.out ./mainB_row 10 >> B.out
	\time -f "\t\ttime: %U seconds" -ao B.out ./mainB_row 50 >> B.out
	\time -f "\t\ttime: %U seconds" -ao B.out ./mainB_row 100 >> B.out
	\time -f "\t\ttime: %U seconds" -ao B.out ./mainB_row 150 >> B.out
	echo "\nMaximum element of row:" >> B.out
	\time -f "\t\ttime: %U seconds" -ao B.out ./mainB_max 10 >> B.out
	\time -f "\t\ttime: %U seconds" -ao B.out ./mainB_max 50 >> B.out
	\time -f "\t\ttime: %U seconds" -ao B.out ./mainB_max 100 >> B.out
	\time -f "\t\ttime: %U seconds" -ao B.out ./mainB_max 150 >> B.out

