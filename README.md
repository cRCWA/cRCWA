# cRCWA

cRCWA is an efficient implementation of the Rigorous Coupled Wave Analysis (RCWA) and Aperiodic Fourier Modal Method (AFMM) techniques. It has been developed for simulating a relatively large set of photonic systems. It is developed at the Centre de Radiofréquences, Optique et Micro-nanoélectronique des Alpes (CROMA) in Grenoble, France and it has been successfully used to study diverse structures, ranging from microstructured solar cells to curved waveguides.

This open source project is distributed under the GPL v.3 license.

## Compile and install

cRCWA is written in C++ and can be compiled on a Unix system using Make and gcc. On a practical standpoint, it has been tested on Linux and on macOS.
You will need the following libraries installed in your system:
- LAPACK, the classic library for linear algebra
- BLAS, the basic linear algebra package, required by LAPACK
- FFTW3, the library to efficiently compute fast Fourier transforms

The source files of cRCWA and the makefile are present in the src directory.
Write down the folders where the libraries are installed in your system, as well as the name of the libraries. For instance, let us suppose that in your system the LAPACK library is stored as `/usr/local/lib/liblapack.so`. In this case, change the configuration section of the `makefile` (in the `src` directory) so that for the LAPACK library you have:

~~~~
LIBLAPACKDIR = /usr/local/lib/
LIBLAPACK = -llapack
~~~~~

When the makefile is correctly configured, you should be able to compile cRCWA by typing `make`, while in the `src` subdirectory.

If the compile is successful, the executables are available in the bin directory.

## Compile cRCWA on macOS

You can use macports to install the LAPACK, BLAS and FFTW3 libraries. The macports tool installs libraries in `/opt/local/lib`. The LAPACK library is installed in a subdirectory of this folder. In this case, the configuration section at the beginning of your makefile may look as follows:

~~~~
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

# ***************************************************************************
~~~~

## User manual

cRCWA is described in an user manual available here: https://github.com/cRCWA/cRCWA/blob/main/manual/manual/afmm.pdf

## Python integration

cRCWA can be used as a Python module...