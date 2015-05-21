CC = gcc
CFLAGS = -Wall -lm -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

# Overall
all: A

clean:
	rm -f *.o *.dat *.out *.pdf mainA mainB mainC

# Object files
ode_rkf45.o: ../7_ode/ode_rkf45.c ../7_ode/ode_rkfe5.h
	$(CC) $(CFLAGS) -c $< -o $@
ode_driver.o: ../7_ode/ode_driver.c ../7_ode/ode_driver.c
	$(CC) $(CFLAGS) -c $< -o $@
ode_integrate.o: ../7_ode/ode_integrate.c ../7_ode/ode_integrate.h
	$(CC) $(CFLAGS) -c $< -o $@
qarc.o: qarc.c qarc.h
qaro.o: qaro.c qaro.h

# Task A
A: A.out

A.out: mainA
	./mainA > A.out
mainA: mainA.c qarc.o qaro.o

