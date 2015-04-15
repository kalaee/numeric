
all: A B
clean: cleanA cleanB
check: checkA checkB

# Project A
A: out.A.txt
	./mainA

out.A.txt: mainA
	./mainA

mainA: mainA.c

checkA: out.A.txt
	@echo "checking target A ..."
	@printf "Hello, World!\n" > correct.txt
	@diff --brief correct.txt out.A.txt
	@echo "target A seems to be ok"
	@rm -f correct.txt

cleanA:
	rm -f mainA out.A.txt

#Project B
B: out.B.txt


out.B.txt: mainB
	./mainB > out.B.txt

mainB: mainB.o hello.o world.o
mainB.o: mainB.c
hello.o: hello.c
world.o: world.c

checkB: out.B.txt
	@echo "checking target B ..."
	@printf "Hello, World!\n" > correct.txt
	@diff --brief correct.txt out.B.txt
	@echo "target B seems to be ok"
	@rm -f correct.txt

cleanB:
	rm -f *.o out.B.txt mainB
