CFLAGS = -Wall -lm -std=gnu99 -g

# All
all: A B C

clean: cleanA cleanB cleanC

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
	
# Project B
C: C.dat plotC.ppl
	pyxplot plotC.ppl

C.dat: mainC
	./mainC

mainC: mainC.c cspline.h

cleanC:
	rm -f C.dat mainC plotC.pdf
