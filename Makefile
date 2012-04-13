CC = g++
CFLAGS = -c -Wall -O3
INCLUDES=
# program name
MAIN = quantwo
# base directory
BASE=$(shell dirname $$(readlink -f Makefile))
# files to be linked to working-directory
FILIN=definitions.tex equation.tex
DIR = src
OBJ0 = main.o finput.o orbital.o matrices.o operators.o kronecker.o term.o globals.o
OBJ = $(patsubst %,$(DIR)/%,$(OBJ0))
SRC = $(OBJ:.o=.cpp)

all : $(MAIN)

$(MAIN) : $(OBJ)
	$(CC) $(OBJ) $(INCLUDES) -o $(MAIN)

clean : 
	rm -rf $(MAIN) $(OBJ) 
veryclean :
	rm -rf $(MAIN) $(OBJ) *.pdf *.aux *.bib *.dvi *.ps *.log *~ $(DIR)/*~  
equation :
	 $(MAIN) input
	 pdflatex equation.tex

base : $(FILIN)
ifneq ($(BASE),$(PWD))
	 @test -e $(MAIN) || ln -s $(BASE)/$(MAIN) .
endif

$(FILIN) :
ifneq ($(BASE),$(PWD))
	 @test -e $@ || ln -s $(BASE)/$@ .
endif

depend: $(SRC)
	 makedepend -- $(INCLUDES) $(CFLAGS) -Y -- $^

# DO NOT DELETE THIS LINE -- make depend needs it

src/main.o: src/utilities.h src/utilities.cpp src/term.h src/product.h
src/main.o: src/product.cpp src/operators.h src/globals.h src/orbital.h
src/main.o: src/matrices.h src/kronecker.h src/sum.h src/sum.cpp src/finput.h
src/finput.o: src/finput.h src/utilities.h src/utilities.cpp src/product.h
src/finput.o: src/product.cpp src/term.h src/operators.h src/globals.h
src/finput.o: src/orbital.h src/matrices.h src/kronecker.h src/sum.h
src/finput.o: src/sum.cpp
src/orbital.o: src/orbital.h src/utilities.h src/utilities.cpp
src/matrices.o: src/matrices.h src/product.h src/utilities.h
src/matrices.o: src/utilities.cpp src/product.cpp src/orbital.h src/globals.h
src/operators.o: src/operators.h src/utilities.h src/utilities.cpp
src/operators.o: src/globals.h src/orbital.h src/matrices.h src/product.h
src/operators.o: src/product.cpp
src/kronecker.o: src/kronecker.h src/orbital.h src/utilities.h
src/kronecker.o: src/utilities.cpp src/globals.h
src/term.o: src/term.h src/utilities.h src/utilities.cpp src/product.h
src/term.o: src/product.cpp src/operators.h src/globals.h src/orbital.h
src/term.o: src/matrices.h src/kronecker.h src/sum.h src/sum.cpp
src/globals.o: src/globals.h
