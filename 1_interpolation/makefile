CFLAGS = -Wall -lm -std=gnu99 -g

# All
all: A B

clean: cleanA cleanB

# Project A
A: A.dat plotA.ppl
	pyxplot plotA.ppl

A.dat: mainA
	./mainA

mainA: mainA.c qspline.h lspline.h

cleanA:
	rm -f A.dat mainA plotA.pdf

# Project B
B: B.dat plotB.ppl
	pyxplot plotB.ppl

B.dat: mainB
	./mainB

mainB: mainB.c qspline.h

cleanB:
	rm -f B.dat mainB plotB.pdf
