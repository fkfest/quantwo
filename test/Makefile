# program name
MAIN = ../quantwo
OBJ = ccsd ip

test : $(OBJ)

$(OBJ) :
	@test -e $@ || $(MAIN) $@.q2
	diff -q $@.tex $@.ref 
	rm -f $@.tex
