CC = gcc
CFLAGS = -Wall -lm -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

# Overall
all: A B C

clean:
	rm -f *.o *.out *.pdf mainA mainB mainC

# Object files
ode_rkf45.o: ../libs/ode_rkf45.c ../libs/ode_rkf45.h
	$(CC) $(CFLAGS) -c $< -o $@
ode_driver.o: ../libs/ode_driver.c ../libs/ode_driver.c
	$(CC) $(CFLAGS) -c $< -o $@
ode_integrate.o: ../libs/ode_integrate.c ../libs/ode_integrate.h
	$(CC) $(CFLAGS) -c $< -o $@
qarc.o: qarc.c qarc.h
qaro.o: qaro.c qaro.h
qasc.o: qasc.c qasc.h
qaso.o: qaso.c qaso.h

# Task A
A: A.out
A.out: mainA
	./mainA > A.out
mainA: mainA.c qarc.o qaro.o

# Task B
B: B.out
B.out: mainB
	./mainB > B.out
mainB: mainB.c qaso.o qasc.o qaro.o

# Task C
C: C.out
C.out: mainC
	./mainC > C.out
mainC: mainC.c qaso.o qasc.o qaro.o qarc.o ode_rkf45.o ode_driver.o ode_integrate.o
