CC = gcc
CFLAGS = -Wall -lm -std=gnu99
CFLAGS += `gsl-config --cflags`
LDLIBS += `gsl-config --libs`

# Overall

all: A B C

clean:
	rm -f *.o *.out *.pdf mainA mainB mainC mainB_givens mainB_qr

# Object files
qr.o: qr.c backsub.h
givens.o: givens.c backsub.h
cholesky.o: cholesky.c cholesky.h
backsub.o: backsub.c backsub.h

# Task A
A: mainA
	./mainA > A.out

mainA: mainA.c qr.o backsub.o

# Task B
B: mainB B.pdf
	./mainB > B.out
mainB: mainB.c givens.o backsub.o
mainB_givens: mainB_givens.c givens.o backsub.o
mainB_qr: mainB_qr.c qr.o backsub.o
time.out: mainB_givens mainB_qr
	echo '# Givens qr' > $@
	for N in `seq 100 10 1000`; do \
	echo "N=$$N"; \time -ao $@ -f "$$N %U" ./mainB_givens $$N; \
	done
	echo '\n\n# Gram-Schmidt qr' >> $@
	for N in `seq 100 10 1000`; do \
	echo "N=$$N"; \time -ao $@ -f "$$N %U" ./mainB_qr $$N; \
	done
B.pdf: B.ppl time.out
	pyxplot B.ppl

# Task C
C: mainC
	./mainC > C.out

mainC: mainC.c cholesky.o backsub.o
