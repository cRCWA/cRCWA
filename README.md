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

## Self tests

cRCWA contains a certain number of self tests that can be run in the `test` directory. They can be useful to check for problems in your install. We encourage to run them often via the `run_all_tests.sh` script.

A failure in a test does not immediately mean that a problem is present, especially when results are compared with respect to a reference. The tests try to check if two results are reasonably close. It is in fact difficult to compare two simulated results for the electric fields. Some small variations may be normal and due to the truncation to a given floating point precision.

## User manual

cRCWA is described in an user manual available here: https://github.com/cRCWA/cRCWA/blob/main/manual/manual/afmm.pdf

## jOptiEx

cRCWA may generate large collections of files in some situations. For instance, when exporting the modes of a multimodal waveguide or structure. JOptiEx is a Java tool that allows to explore very rapidly a collection of files and represent them rapidly. It is not meant as a full-fledged scientific representation tool, but it is more a tool to explore the results to select those to be represented more carefully.

jOptiEx is a Java program. To run it, you need to have a JRE available in your computer and type `java -jar jOptiEx.jar` from the directory where `jOptiEx.jar` is present.


## Python integration

cRCWA can be used as a Python module...