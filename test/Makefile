# program name
MAIN = ../quantwo
OBJ = ccsd ugccsd ugdcsd sums ip mr names genwick connections

test : $(OBJ)
			 diff -q ugccsd.alg ugccsd_alg.ref
			 diff -q ugdcsd.jl ugdcsd_jl.ref
			 rm ugdcsd.jl ugccsd.alg

$(OBJ) :
	@test -e $@ || $(MAIN) $@.q2
	diff -q $@.tex $@.ref 
	rm -f $@.tex
