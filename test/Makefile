# program name
MAIN = ../quantwo
OBJ = ccsd ugccsd ugdcsd sums ip mr names genwick connections

test : $(OBJ)

$(OBJ) :
	@test -e $@ || $(MAIN) $@.q2
	diff -q $@.tex $@.ref 
	rm -f $@.tex
