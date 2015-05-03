CC = gcc
CFLAGS = -Wall -lm -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

# Overall

all: A

clean:
	rm -f *.o *.out mainA mainB_cyc mainB_row mainB_max

# Object files
jacobi.o: jacobi.c

# Task A
A: mainA
	./mainA > A.out

mainA: mainA.c jacobi.o jacobi.h
mainB_cyc: mainB_cyc.c jacobi.o jacobi.h
mainB_row: mainB_row.c jacobi.o jacobi.h
mainB_max: mainB_max.c jacobi.o jacobi.h
mainB: mainB.c jacobi.c jacobi.h

B: mainB_cyc mainB_row mainB_max
	\time -f "%U" -o time_cyc.out ./mainB_cyc 10 > B_cyc.out
	\time -f "%U" -ao time_cyc.out ./mainB_cyc 20 >> B_cyc.out
	\time -f "%U" -ao time_cyc.out ./mainB_cyc 40 >> B_cyc.out
	\time -f "%U" -ao time_cyc.out ./mainB_cyc 60 >> B_cyc.out
	\time -f "%U" -ao time_cyc.out ./mainB_cyc 80 >> B_cyc.out
	\time -f "%U" -ao time_cyc.out ./mainB_cyc 100 >> B_cyc.out
	\time -f "%U" -ao time_cyc.out ./mainB_cyc 120 >> B_cyc.out
	\time -f "%U" -ao time_cyc.out ./mainB_cyc 140 >> B_cyc.out
	\time -f "%U" -ao time_cyc.out ./mainB_cyc 160 >> B_cyc.out
	\time -f "%U" -ao time_cyc.out ./mainB_cyc 180 >> B_cyc.out
	\time -f "%U" -ao time_cyc.out ./mainB_cyc 200 >> B_cyc.out
	\time -f "%U" -ao time_cyc.out ./mainB_cyc 220 >> B_cyc.out
	\time -f "%U" -ao time_cyc.out ./mainB_cyc 240 >> B_cyc.out
	\time -f "%U" -ao time_cyc.out ./mainB_cyc 260 >> B_cyc.out
	\time -f "%U" -ao time_cyc.out ./mainB_cyc 280 >> B_cyc.out
	\time -f "%U" -ao time_cyc.out ./mainB_cyc 300 >> B_cyc.out
	\time -f "%U" -ao time_cyc.out ./mainB_cyc 320 >> B_cyc.out
	\time -f "%U" -ao time_cyc.out ./mainB_cyc 340 >> B_cyc.out
	\time -f "%U" -ao time_cyc.out ./mainB_cyc 360 >> B_cyc.out
	\time -f "%U" -ao time_cyc.out ./mainB_cyc 400 >> B_cyc.out
	\time -f "%U" -ao time_cyc.out ./mainB_cyc 440 >> B_cyc.out
	\time -f "%U" -ao time_cyc.out ./mainB_cyc 480 >> B_cyc.out
	\time -f "%U" -o time_row.out ./mainB_row 10 > B_row.out
	\time -f "%U" -ao time_row.out ./mainB_row 20 >> B_row.out
	\time -f "%U" -ao time_row.out ./mainB_row 40 >> B_row.out
	\time -f "%U" -ao time_row.out ./mainB_row 60 >> B_row.out
	\time -f "%U" -ao time_row.out ./mainB_row 80 >> B_row.out
	\time -f "%U" -ao time_row.out ./mainB_row 100 >> B_row.out
	\time -f "%U" -ao time_row.out ./mainB_row 120 >> B_row.out
	\time -f "%U" -ao time_row.out ./mainB_row 140 >> B_row.out
	\time -f "%U" -ao time_row.out ./mainB_row 160 >> B_row.out
	\time -f "%U" -ao time_row.out ./mainB_row 180 >> B_row.out
	\time -f "%U" -ao time_row.out ./mainB_row 200 >> B_row.out
	\time -f "%U" -ao time_row.out ./mainB_row 220 >> B_row.out
	\time -f "%U" -ao time_row.out ./mainB_row 240 >> B_row.out
	\time -f "%U" -ao time_row.out ./mainB_row 260 >> B_row.out
	\time -f "%U" -ao time_row.out ./mainB_row 280 >> B_row.out
	\time -f "%U" -ao time_row.out ./mainB_row 300 >> B_row.out
	\time -f "%U" -ao time_row.out ./mainB_row 320 >> B_row.out
	\time -f "%U" -ao time_row.out ./mainB_row 340 >> B_row.out
	\time -f "%U" -ao time_row.out ./mainB_row 360 >> B_row.out
	\time -f "%U" -ao time_row.out ./mainB_row 400 >> B_row.out
	\time -f "%U" -ao time_row.out ./mainB_row 440 >> B_row.out
	\time -f "%U" -ao time_row.out ./mainB_row 480 >> B_row.out
	\time -f "%U" -o time_max.out ./mainB_max 10 > B_max.out
	\time -f "%U" -ao time_max.out ./mainB_max 20 >> B_max.out
	\time -f "%U" -ao time_max.out ./mainB_max 40 >> B_max.out
	\time -f "%U" -ao time_max.out ./mainB_max 60 >> B_max.out
	\time -f "%U" -ao time_max.out ./mainB_max 80 >> B_max.out
	\time -f "%U" -ao time_max.out ./mainB_max 100 >> B_max.out
	\time -f "%U" -ao time_max.out ./mainB_max 120 >> B_max.out
	\time -f "%U" -ao time_max.out ./mainB_max 140 >> B_max.out
	\time -f "%U" -ao time_max.out ./mainB_max 160 >> B_max.out
	\time -f "%U" -ao time_max.out ./mainB_max 180 >> B_max.out
	\time -f "%U" -ao time_max.out ./mainB_max 200 >> B_max.out
	\time -f "%U" -ao time_max.out ./mainB_max 220 >> B_max.out
	\time -f "%U" -ao time_max.out ./mainB_max 240 >> B_max.out
	\time -f "%U" -ao time_max.out ./mainB_max 260 >> B_max.out
	\time -f "%U" -ao time_max.out ./mainB_max 280 >> B_max.out
	\time -f "%U" -ao time_max.out ./mainB_max 300 >> B_max.out
	\time -f "%U" -ao time_max.out ./mainB_max 320 >> B_max.out
	\time -f "%U" -ao time_max.out ./mainB_max 340 >> B_max.out
	\time -f "%U" -ao time_max.out ./mainB_max 360 >> B_max.out
	\time -f "%U" -ao time_max.out ./mainB_max 400 >> B_max.out
	\time -f "%U" -ao time_max.out ./mainB_max 440 >> B_max.out
	\time -f "%U" -ao time_max.out ./mainB_max 480 >> B_max.out

