# ***************************************************************************
# * GNU makefile for compiling the program on a UNIX system                 *
# *                                                                         *
# * This program needs Lapack, Blas as well as f2c and fftw3                *
# ***************************************************************************
# 
# This file is part of cRCWA.
# 
#   cRCWA is free software: you can redistribute it and/or modify it under the
#   terms of the GNU General Public License as published by the Free Software
#   Foundation, either version 3 of the License, or (at your option) any later
#   version.
#   
#   cRCWA is distributed in the hope that it will be useful, but WITHOUT ANY
#   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
#   FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
#   details.
#   
#   You should have received a copy of the GNU General Public License along
#   with cRCWA. If not, see <https://www.gnu.org/licenses/>. 
# 
#   Davide Bucci, 2008-2025
#   Jérôme Michallon, 2012-2014
# */


# ***************************************************************************
# CONFIGURATION SECTION
# ***************************************************************************

# This is the directory where the library files are provided. It is useful
# if the LAPACK, BLAS and FFTW3 libraries are stored in the same place in
# your system. This is very often the case.

LIBDIR =  /opt/local/lib

# This is where the LAPACK lib can be found
LIBLAPACKDIR = $(LIBDIR)/lapack
# The option to include the LAPACK library. You rarely need to change it, as
# the name of the library is very often just "lapack".
LIBLAPACK = -llapack

# This is where the BLAS lib can be found
LIBBLASDIR = $(LIBDIR)/lapack
# The option to include the BLAS library. You will have to change this if you
# are using a library having a different name. For example, a common choice
# with ATLAS may be LIBBLAS = -latlas.
LIBBLAS = -lblas


# This is where the fftw3 lib can be found
LIBFFTW3DIR = $(LIBDIR)
# The option to include the FFTW3 library:
LIBFFTW3 = -lfftw3

# This is where the Python include file Python.h is present:
PYTHONINCDIR = /System/Volumes/Data/opt/local/Library/Frameworks/Python.framework/Versions/3.13/include/python3.13/

# This is the Python library to be included. In macOS the file can be 
# called like libpython3.15.dylib (for a dynamically linked library).
# If you use, let's say, Python 3.11, you should change the version
# accordingly. cRCWA requires Python 3 in any case.

PYTHONLIB = -L/opt/local/Library/Frameworks/Python.framework/Versions/3.15/lib/ -lpython3.15

# ***************************************************************************

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

# Put pycrcwa after o2f if you want the Python library to be automatically
# built with make all.
all: setrev crcwa  slice o2f
	mv crcwa ../bin
	mv slice ../bin
	mv o2g ../bin

# monocore
crcwa: $(OFILES) $(HEADERS)
	$(LD) $(OFILES) -L$(LIBLAPACKDIR) $(LIBLAPACK) -L$(LIBBLASDIR) $(LIBBLAS) -L$(LIBFFTW3DIR) -lfftw3 $(LIBTHREAD)   -o crcwa

pycrcwa:
	$(CXX) -shared -L$(LIBLAPACKDIR) $(LIBLAPACK) -L$(LIBBLASDIR) $(LIBBLAS) -L$(LIBFFTW3DIR) -lfftw3 $(LIBTHREAD) $(PYTHONLIB) -o pycRCWA.so \
	 $(OFILES) -I$(PYTHONINCDIR) afmm_python_module.cpp -fPIC
	 mv pycRCWA.so ../python



slice: $(SLOFILES)
	$(LD) $(SLOFILES) -o slice

o2f: $(O2GFILES)
	$(LD) $(O2GFILES) -o o2g

clean:
	rm *.o
