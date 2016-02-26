CC = g++ 
PROFILE =
#PROFILE = -pg
#PROFILE = -g
CFLAGS := -c -Wall -Wextra -pedantic -std=c++11 -Ofast $(PROFILE)
LDFLAGS = $(PROFILE)
#comment out to deactivate debug and asserts
#CFLAGS := $(CFLAGS) -D NDEBUG
#use rational numbers from boost
#CFLAGS := $(CFLAGS) -D _RATIONAL
INCLUDES=
# program name
MAIN = quantwo
# OS type
UNAME_S := $(shell uname -s)
# base directory
BASE=$(PWD)
ifeq ($(UNAME_S),Linux)
  BASE=$(shell dirname $$(readlink -f Makefile))
endif
ifeq ($(findstring CYGWIN, $(UNAME_S)),CYGWIN)
  BASE=$(shell dirname $$(readlink -f Makefile))
endif
ifeq ($(UNAME_S),Darwin)
  # uses greadlink from coreutils
  BASE=$(shell dirname $$(greadlink -f Makefile))
endif
in=input
out=equation
# files to be linked to working-directory
FILIN=definitions.tex $(out).tex
DIR = src
OBJ0 = main.o tensor.o action.o expression.o factorizer.o unigraph.o inpline.o finput.o equation.o lexic.o work.o orbital.o matrix.o operators.o kronecker.o term.o utilities.o globals.o
OBJ = $(patsubst %,$(DIR)/%,$(OBJ0))
SRC = $(OBJ:.o=.cpp)

all : $(MAIN)

$(MAIN) : $(OBJ)
	$(CC) $(LDFLAGS) $(OBJ) $(INCLUDES) -o $(MAIN)

%.o : %.cpp 
	$(CC) $(CFLAGS) $< -o $@

clean : 
	rm -rf $(MAIN) $(OBJ) 
veryclean :
	git clean -dfx  
equation :
	 $(MAIN) $(in).q2
	 @test -e $(out).tex || cp equation.tex $(out).tex
	 pdflatex "\newcommand\QuantwoInputFileName{$(in)}\input{$(out)}"
	 rm $(out).log $(out).aux

base : $(FILIN)
ifneq ($(BASE),$(PWD))
	 @test -e $(MAIN) || ln -s $(BASE)/$(MAIN) .
endif

$(FILIN) :
ifneq ($(BASE),$(PWD))
	 @test -e $@ || ln -s $(BASE)/$@ .
endif

depend: 
	 makedepend -- $(INCLUDES) $(CFLAGS) -Y -- $(SRC) 

# DO NOT DELETE THIS LINE -- make depend needs it

src/main.o: src/utilities.h src/globals.h src/term.h src/types.h
src/main.o: src/product.h src/product.cpp src/arrays.h src/arrays.cpp
src/main.o: src/operators.h src/orbital.h src/inpline.h src/matrix.h
src/main.o: src/sum.h src/sum.cpp src/kronecker.h src/evertices.h
src/main.o: src/finput.h src/equation.h src/lexic.h src/work.h src/unigraph.h
src/main.o: src/factorizer.h src/tensor.h src/action.h src/expression.h
src/tensor.o: src/tensor.h src/globals.h src/utilities.h src/types.h
src/tensor.o: src/product.h src/product.cpp src/arrays.h src/arrays.cpp
src/action.o: src/action.h src/globals.h src/utilities.h src/types.h
src/action.o: src/product.h src/product.cpp src/arrays.h src/arrays.cpp
src/action.o: src/tensor.h
src/expression.o: src/expression.h src/globals.h src/utilities.h src/types.h
src/expression.o: src/product.h src/product.cpp src/arrays.h src/arrays.cpp
src/expression.o: src/tensor.h src/action.h src/sum.h src/sum.cpp
src/factorizer.o: src/factorizer.h src/globals.h src/types.h src/product.h
src/factorizer.o: src/utilities.h src/product.cpp src/arrays.h src/arrays.cpp
src/factorizer.o: src/orbital.h src/inpline.h src/sum.h src/sum.cpp
src/factorizer.o: src/term.h src/operators.h src/matrix.h src/kronecker.h
src/factorizer.o: src/evertices.h src/tensor.h src/action.h src/expression.h
src/unigraph.o: src/unigraph.h src/types.h src/product.h src/utilities.h
src/unigraph.o: src/globals.h src/product.cpp src/arrays.h src/arrays.cpp
src/unigraph.o: src/term.h src/operators.h src/orbital.h src/inpline.h
src/unigraph.o: src/matrix.h src/sum.h src/sum.cpp src/kronecker.h
src/unigraph.o: src/evertices.h
src/inpline.o: src/inpline.h src/utilities.h src/globals.h
src/finput.o: src/finput.h src/utilities.h src/globals.h src/product.h
src/finput.o: src/product.cpp src/term.h src/types.h src/arrays.h
src/finput.o: src/arrays.cpp src/operators.h src/orbital.h src/inpline.h
src/finput.o: src/matrix.h src/sum.h src/sum.cpp src/kronecker.h
src/finput.o: src/evertices.h src/equation.h src/lexic.h
src/equation.o: src/equation.h src/utilities.h src/globals.h src/product.h
src/equation.o: src/product.cpp src/term.h src/types.h src/arrays.h
src/equation.o: src/arrays.cpp src/operators.h src/orbital.h src/inpline.h
src/equation.o: src/matrix.h src/sum.h src/sum.cpp src/kronecker.h
src/equation.o: src/evertices.h src/lexic.h
src/lexic.o: src/lexic.h src/utilities.h src/globals.h src/product.h
src/lexic.o: src/product.cpp src/types.h src/arrays.h src/arrays.cpp
src/work.o: src/work.h src/utilities.h src/globals.h src/product.h
src/work.o: src/product.cpp src/operators.h src/types.h src/arrays.h
src/work.o: src/arrays.cpp src/orbital.h src/inpline.h src/matrix.h src/sum.h
src/work.o: src/sum.cpp src/kronecker.h src/evertices.h src/term.h
src/work.o: src/unigraph.h src/finput.h src/equation.h src/lexic.h
src/work.o: src/factorizer.h src/tensor.h src/action.h src/expression.h
src/orbital.o: src/orbital.h src/utilities.h src/globals.h src/product.h
src/orbital.o: src/product.cpp src/inpline.h
src/matrix.o: src/matrix.h src/globals.h src/types.h src/product.h
src/matrix.o: src/utilities.h src/product.cpp src/arrays.h src/arrays.cpp
src/matrix.o: src/orbital.h src/inpline.h src/sum.h src/sum.cpp
src/matrix.o: src/kronecker.h src/evertices.h
src/operators.o: src/operators.h src/utilities.h src/globals.h src/types.h
src/operators.o: src/product.h src/product.cpp src/arrays.h src/arrays.cpp
src/operators.o: src/orbital.h src/inpline.h src/matrix.h src/sum.h
src/operators.o: src/sum.cpp src/kronecker.h src/evertices.h src/term.h
src/kronecker.o: src/kronecker.h src/orbital.h src/utilities.h src/globals.h
src/kronecker.o: src/product.h src/product.cpp src/inpline.h
src/term.o: src/term.h src/utilities.h src/globals.h src/types.h
src/term.o: src/product.h src/product.cpp src/arrays.h src/arrays.cpp
src/term.o: src/operators.h src/orbital.h src/inpline.h src/matrix.h
src/term.o: src/sum.h src/sum.cpp src/kronecker.h src/evertices.h
src/utilities.o: src/utilities.h src/globals.h
src/globals.o: src/globals.h
