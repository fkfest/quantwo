# program name
MAIN = ../quantwo
OBJ = ccsd ugccsd ugdcsd dcsd dcsd-2 uccsd uccsdt udc-ccsdt dc-ccsdt tc-ccsd tc-uccsd sums ip mr names genwick connections

test : $(OBJ)
			 @echo "All tests passed."

$(OBJ) :
	@test -e $@ || $(MAIN) $@.q2
	diff -q $@.tex $@.ref 
	rm -f $@.tex
	@test ! -e $@.alg || (diff -q $@.alg $@_alg.ref && rm $@.alg)
	@test ! -e $@.jl || (diff -q $@.jl $@_jl.ref && rm $@.jl)
