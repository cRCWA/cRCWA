# ***************************************************************************
# * GNU makefile for compiling the program on a UNIX system                 *
# *                                                                         *
# * This program needs Lapack, Blas as well as f2c and fftw3                *
# ***************************************************************************

# I have chosen to put all the useful files in the single 
# directory specified by the $(LIBDIR) variable.

# This is the directory where the library files are provided:
LIBDIR =  /opt/local/lib

# This is where the LAPACK lib can be found
LIBLAPACK = $(LIBDIR)/lapack

# This is where the BLAS lib can be found
LIBBLAS = $(LIBDIR)/lapack

# This is where the fftw3 lib can be found
LIBFFTW3 = $(LIBDIR)

# For multi-thread operation, we need:
LIBTHREAD = -lpthread

# Customize the options following the name conventions in your library.
# If in doubt, check with the nm command on the library file installed in your
# system. Variants may require a leading or final underscore. Remember that C
# often adds automatically an underscore at the beginning of the function name.
# Therefore if you see a name like _zgtri_ you should use 
# LAPACK_ADD_FINAL_UNDERSCORE
#
# The available options are:
# LAPACK_ADD_FINAL_UNDERSCORE
# LAPACK_ADD_LEADING_UNDERSCORE
# LAPACK_ADD_BOTH_UNDERSCORES
#
# BLAS_ADD_FINAL_UNDERSCORE
# BLAS_ADD_LEADING_UNDERSCORE
# BLAS_ADD_BOTH_UNDERSCORES

OPTIONS = -D LAPACK_ADD_FINAL_UNDERSCORE -D BLAS_ADD_FINAL_UNDERSCORE



.SUFFIXES:
.SUFFIXES: .cpp .o
CC = gcc
LD = g++
SHELL = /bin/sh

# These are the flags passed to the C++ compiler. Change accordingly to your
# system
CFLAGS = -g -O3
# CFLAGS = -pedantic -Wall --unroll-all-loops -O3  

.PHONY : setrev

HEADERS =  block_matrix.h  fortran_types.h parsefile.h structure.h section.h commands.h finterface.h parsercl.h matrixDev.h matrixdevsym.h parallelize.h
OFILES =   block_matrix.o main.o matrixDev.o parsefile.o structure.o section.o commands_structure.o commands_section.o commands_propagation.o commands_other.o finterface.o parsercl.o matrixdevsym.o parallelize.o toeplitz.o
SLOFILES = slice.o
O2GFILES = optiwave2gnuplot.o

# compile rule
.cpp.o:
	$(CC) $(OPTIONS) -c $(CFLAGS) $*.cpp 

all: setrev afmm  slice o2f 
	mv afmm ../bin
	mv slice ../bin
	mv o2g ../bin

# monocore
afmm: $(OFILES) $(HEADERS)
	$(LD) $(OFILES) -L$(LIBLAPACK) -llapack -L$(LIBBLAS) -lblas -L$(LIBFFTW3) -lfftw3 $(LIBTHREAD)   -o afmm

slice: $(SLOFILES)
	$(LD) $(SLOFILES) -o slice

o2f: $(O2GFILES)
	$(LD) $(O2GFILES) -o o2g

clean:
	rm *.o
