CC = g++ 
PROFILE =
#PROFILE = -pg
CFLAGS := -c -Wall -O2 $(PROFILE)
LDFLAGS = $(PROFILE)
#use rational numbers from boost
CFLAGS := $(CFLAGS) -D _RATIONAL
INCLUDES=
# program name
MAIN = quantwo
# base directory
BASE=$(shell dirname $$(readlink -f Makefile))
# files to be linked to working-directory
FILIN=definitions.tex equation.tex
DIR = src
OBJ0 = main.o finput.o work.o orbital.o matrices.o operators.o kronecker.o term.o utilities.o globals.o
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
	rm -rf $(MAIN) $(OBJ) *.pdf *.aux *.bib *.dvi *.ps *.log *~ $(DIR)/*~  
equation :
	 $(MAIN) input.q2
	 pdflatex equation.tex

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

src/main.o: src/utilities.h src/globals.h src/term.h src/product.h
src/main.o: src/product.cpp src/operators.h src/orbital.h src/matrices.h
src/main.o: src/kronecker.h src/sum.h src/sum.cpp src/finput.h src/work.h
src/finput.o: src/finput.h src/utilities.h src/globals.h src/product.h
src/finput.o: src/product.cpp src/term.h src/operators.h src/orbital.h
src/finput.o: src/matrices.h src/kronecker.h src/sum.h src/sum.cpp
src/work.o: src/work.h src/utilities.h src/globals.h src/product.h
src/work.o: src/product.cpp src/operators.h src/orbital.h src/matrices.h
src/work.o: src/kronecker.h src/sum.h src/sum.cpp src/term.h
src/orbital.o: src/orbital.h src/utilities.h src/globals.h
src/matrices.o: src/matrices.h src/product.h src/utilities.h src/globals.h
src/matrices.o: src/product.cpp src/orbital.h
src/operators.o: src/operators.h src/utilities.h src/globals.h src/orbital.h
src/operators.o: src/matrices.h src/product.h src/product.cpp
src/kronecker.o: src/kronecker.h src/orbital.h src/utilities.h src/globals.h
src/term.o: src/term.h src/utilities.h src/globals.h src/product.h
src/term.o: src/product.cpp src/operators.h src/orbital.h src/matrices.h
src/term.o: src/kronecker.h src/sum.h src/sum.cpp
src/utilities.o: src/utilities.h src/globals.h
src/globals.o: src/globals.h
