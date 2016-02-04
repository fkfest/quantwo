CC = g++ 
PROFILE =
#PROFILE = -pg
CFLAGS := -c -Wall -Wextra -pedantic -Ofast $(PROFILE)
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
OBJ0 = main.o tensor.o action.o expression.o factorizer.o inpline.o finput.o equation.o lexic.o work.o orbital.o matrices.o operators.o kronecker.o term.o utilities.o globals.o
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
src/main.o: src/product.h src/product.cpp src/operators.h src/orbital.h
src/main.o: src/matrices.h src/inpline.h src/sum.h src/sum.cpp
src/main.o: src/kronecker.h src/finput.h src/equation.h src/lexic.h
src/main.o: src/work.h src/factorizer.h src/tensor.h src/arrays.h
src/main.o: src/arrays.cpp src/action.h src/expression.h
src/tensor.o: src/tensor.h src/globals.h src/utilities.h src/types.h
src/tensor.o: src/product.h src/product.cpp src/arrays.h src/arrays.cpp
src/action.o: src/action.h src/globals.h src/utilities.h src/types.h
src/action.o: src/product.h src/product.cpp src/tensor.h src/arrays.h
src/action.o: src/arrays.cpp
src/expression.o: src/expression.h src/globals.h src/utilities.h src/types.h
src/expression.o: src/product.h src/product.cpp src/tensor.h src/arrays.h
src/expression.o: src/arrays.cpp src/action.h src/sum.h src/sum.cpp
src/factorizer.o: src/factorizer.h src/globals.h src/types.h src/product.h
src/factorizer.o: src/utilities.h src/product.cpp src/orbital.h src/inpline.h
src/factorizer.o: src/sum.h src/sum.cpp src/term.h src/operators.h
src/factorizer.o: src/matrices.h src/kronecker.h src/tensor.h src/arrays.h
src/factorizer.o: src/arrays.cpp src/action.h src/expression.h
src/inpline.o: src/inpline.h src/utilities.h src/globals.h
src/finput.o: src/finput.h src/utilities.h src/globals.h src/product.h
src/finput.o: src/product.cpp src/term.h src/types.h src/operators.h
src/finput.o: src/orbital.h src/matrices.h src/inpline.h src/sum.h
src/finput.o: src/sum.cpp src/kronecker.h src/equation.h src/lexic.h
src/equation.o: src/equation.h src/utilities.h src/globals.h src/product.h
src/equation.o: src/product.cpp src/term.h src/types.h src/operators.h
src/equation.o: src/orbital.h src/matrices.h src/inpline.h src/sum.h
src/equation.o: src/sum.cpp src/kronecker.h src/lexic.h
src/lexic.o: src/lexic.h src/utilities.h src/globals.h src/product.h
src/lexic.o: src/product.cpp src/types.h
src/work.o: src/work.h src/utilities.h src/globals.h src/product.h
src/work.o: src/product.cpp src/operators.h src/types.h src/orbital.h
src/work.o: src/matrices.h src/inpline.h src/sum.h src/sum.cpp
src/work.o: src/kronecker.h src/term.h src/finput.h src/equation.h
src/work.o: src/lexic.h src/factorizer.h src/tensor.h src/arrays.h
src/work.o: src/arrays.cpp src/action.h src/expression.h
src/orbital.o: src/orbital.h src/utilities.h src/globals.h src/product.h
src/orbital.o: src/product.cpp
src/matrices.o: src/matrices.h src/globals.h src/types.h src/product.h
src/matrices.o: src/utilities.h src/product.cpp src/orbital.h src/inpline.h
src/matrices.o: src/sum.h src/sum.cpp src/kronecker.h
src/operators.o: src/operators.h src/utilities.h src/globals.h src/types.h
src/operators.o: src/product.h src/product.cpp src/orbital.h src/matrices.h
src/operators.o: src/inpline.h src/sum.h src/sum.cpp src/kronecker.h
src/operators.o: src/term.h
src/kronecker.o: src/kronecker.h src/orbital.h src/utilities.h src/globals.h
src/kronecker.o: src/product.h src/product.cpp
src/term.o: src/term.h src/utilities.h src/globals.h src/types.h
src/term.o: src/product.h src/product.cpp src/operators.h src/orbital.h
src/term.o: src/matrices.h src/inpline.h src/sum.h src/sum.cpp
src/term.o: src/kronecker.h
src/utilities.o: src/utilities.h src/globals.h
src/globals.o: src/globals.h
