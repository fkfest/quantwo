# program name
MAIN = ../quantwo
OBJ = ccsd ugccsd ugdcsd uccsd uccsdt udc-ccsdt tc-ccsd tc-uccsd sums ip mr names genwick connections

test : $(OBJ)
			 diff -q ugccsd.alg ugccsd_alg.ref
			 diff -q ugdcsd.jl ugdcsd_jl.ref
			 diff -q uccsd.jl uccsd_jl.ref
			 diff -q uccsdt.jl uccsdt_jl.ref
			 diff -q udc-ccsdt.jl udc-ccsdt_jl.ref
			 diff -q tc-ccsd.alg tc-ccsd_alg.ref
			 diff -q tc-uccsd.alg tc-uccsd_alg.ref
			 rm ugdcsd.jl ugccsd.alg uccsd.jl uccsdt.jl udc-ccsdt.jl tc-ccsd.alg tc-uccsd.alg

$(OBJ) :
	@test -e $@ || $(MAIN) $@.q2
	diff -q $@.tex $@.ref 
	rm -f $@.tex
