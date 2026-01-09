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

include ../configure.inc

.SUFFIXES:
.SUFFIXES: .cpp .o


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
