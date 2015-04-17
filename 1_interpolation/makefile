CFLAGS = -Wall -lm -std=gnu99 -g

A: A.dat plotA.ppl
	pyxplot plotA.ppl

A.dat: splineA
	./splineA

splineA: splineA.c qspline.h lspline.h

clean:
	rm -f *.dat splineA *.o *.pdf
