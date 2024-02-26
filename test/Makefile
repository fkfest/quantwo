# program name
MAIN = ../quantwo
OBJ = ccsd ugccsd ugdcsd uccsd uccsdt udc-ccsdt tc-ccsd tc-uccsd sums ip mr names genwick connections

test : $(OBJ)

$(OBJ) :
	@test -e $@ || $(MAIN) $@.q2
	diff -q $@.tex $@.ref 
	rm -f $@.tex
	@test ! -e $@.alg || (diff -q $@.alg $@_alg.ref && rm $@.alg)
	@test ! -e $@.jl || (diff -q $@.jl $@_jl.ref && rm $@.jl)
