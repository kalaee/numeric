CFLAGS = -Wall -lm -std=gnu99 -g

A: A.dat plotA.ppl
	pyxplot plotA.ppl

A.dat: splineA
	./splineA

splineA: splineA.c lspline.o qspline_eval.o qspline_alloc.o qspline_free.o

lspline.o: lspline.c
qspline_eval.o: qspline_eval.c
qspline_alloc.o: qspline_alloc.c
qspline_free.o: qspline_free.c

clean:
	rm -f *.dat splineA *.o
