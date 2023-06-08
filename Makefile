.SUFFIXES : .o .cpp
# compiler and flags
CC     = icc
LINK   = $(CC) -Bstatic
CFLAGS = -O3 $(DEBUG) #-Wno-unused-result

#
OFLAGS = -O3 $(DEBUG)
INC    = $(FFTINC) $(LPKINC) $(USRINC) $(SPGINC) $(GSLINC)
LIB    = $(FFTLIB) $(LPKLIB) $(USRLIB) $(SPGLIB) $(GSLLIB)

# fftw 3 library
FFTINC    = -I/opt/software/fftw/include
FFTLIB    = -L/opt/software/fftw/lib -lfftw3

# Lapack library
#LPKINC = -I/opt/clapack/3.2.1/include
#LPKLIB = -L/opt/clapack/3.2.1/lib -lclapack -lblas -lf2c -lm

# spglib, used to get the irreducible q-points
#SPGINC = -I/opt/spglib/0.7.1/include
#SPGLIB = -L/opt/spglib/0.7.1/lib -lsymspg

# gsl ib, used to get the irreducible q-points
#GSLINC = -I/opt/libs/gsl/include
#GSLLIB = -L/opt/libs/gsl/lib -lgsl -lgslcblas

# Debug flags
#DEBUG = -O -g -DDEBUG
#====================================================================
# executable name
ROOT   = vacf
EXE    = $(ROOT)
#====================================================================
# source and rules
SRC = $(wildcard *.cpp *.c)
OBJ = $(SRC:.cpp=.o)

#====================================================================
all: ver ${EXE}

${EXE}: $(OBJ)
	$(LINK) $(OFLAGS) $(OBJ) $(LIB) -o $@

clean: 
	rm -f *.o *~ *.mod ${EXE}

tar:
	rm -f ${ROOT}.tar; tar -czvf ${ROOT}.tar.gz *.cpp  *.h Makefile README*

ver:
	@echo "#define VERSION `git log|grep '^commit'|wc -l`" > version.h

.f.o:
	$(FC) $(FFLAGS) $(FREE) $(MPI) ${INC} -c $<
.f90.o:
	$(FC) $(FFLAGS) $(FREE) $(MPI) ${INC} -c $<
.c.o:
	$(CC) $(CFLAGS) -c $<
.cpp.o:
	$(CC) $(CFLAGS) $(INC) -c $<
