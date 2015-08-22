# program name
MAIN = ../quantwo
OBJ = ccsd sums ip mr names genwick

test : $(OBJ)

$(OBJ) :
	@test -e $@ || $(MAIN) $@.q2
	diff -q $@.tex $@.ref 
	rm -f $@.tex
