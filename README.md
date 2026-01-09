# cRCWA

cRCWA is an efficient implementation of the Rigorous Coupled Wave Analysis (RCWA) and Aperiodic Fourier Modal Method (AFMM) techniques. It has been developed for simulating a relatively large set of photonic systems. It is developed at the Centre de Radiofréquences, Optique et Micro-nanoélectronique des Alpes (CROMA) in Grenoble, France and it has been successfully used to study diverse structures, ranging from microstructured solar cells to curved waveguides.

This open source project is distributed under the GPL v.3 license.

## Compile and install on a Unix system

cRCWA is written in C++ and can be compiled on a Unix system using GNU Make and gcc. On a practical standpoint, it has been tested on Linux and on macOS.
You will need the following libraries installed in your system:
- LAPACK, the classic library for linear algebra
- BLAS, the basic linear algebra package, required by LAPACK
- FFTW3, the library to efficiently compute fast Fourier transforms

Write down the folders where the libraries are installed in your system, as well as the name of the libraries. For instance, let us suppose that in your system the LAPACK library is stored as `/usr/local/lib/liblapack.so` (the extension `.dylib` is sometimes used on macOS) . In this case, change the configuration section in the `configure.inc` file, so that you have:

~~~~
LIBLAPACKDIR = /usr/local/lib/
LIBLAPACK = -llapack
~~~~~

When the makefile is correctly configured, you should be able to compile cRCWA by typing `make`, while in the project main subdirectory:

~~~~
make
~~~~

If the compile is successful, the executables are available in the `bin` directory. The Python library is put in the `python` directory.

The performances of cRCWA strongly depend on the performances of the BLAS library installed in your system. Therefore, if you want to obtain the fastest calculations possible, make sure you have an optimized version of BLAS in your system. LAPACK ships with a base implementation of BLAS will work fine for test purposes, but will not deliver the best possible performance.  It is not uncommon to cut calculation times by a factor of 2 or 3 shifting from a basic BLAS implementation to an optimized one. Popular and convenient choices may be OpenBLAS or ATLAS.

Some implementations of BLAS or LAPACK may be configured to add an underscore tot the function name. For instance, one of the most important LAPACK functions for cRCWA is `zgeev`. It calculates the eigenvalues and eigenvectors of a non-symmetric matrix of complex elements and is used to calculate the propagation modes. Compilers usually add a leading underscore to the function names. In some cases, another underscore is added at the end of the name. The function `zgeev` therefore appears as `_zgeev_`. For instance, running the `nm` command on macOS we get:
~~~~
% nm /opt/local/lib/lapack/liblapack.dylib |grep zgeev
00000000003146c0 T _zgeev_
000000000071bff0 T _zgeev_64_
0000000000315328 T _zgeevx_
000000000071cc4c T _zgeevx_64_
~~~~
The following option should be provided to add the final underscore every time a LAPACK or BLAS function is called in the code:

~~~~
OPTIONS = -D LAPACK_ADD_FINAL_UNDERSCORE -D BLAS_ADD_FINAL_UNDERSCORE
~~~~

The options `LAPACK_ADD_LEADING_UNDERSCORE` and `LAPACK_ADD_BOTH_UNDERSCORES` are useful in the very rare cases where the names exported by the libraries contain a leading double underscore.

## Compile cRCWA on macOS

You can use macports to install the LAPACK, BLAS and FFTW3 libraries. The macports tool installs libraries in `/opt/local/lib`. The LAPACK library is installed in a subdirectory of this folder.

## Self tests

cRCWA contains a certain number of self tests that can be run in the `test` directory. They can be useful to check for problems in your install. We encourage to run them often via the `run_all_tests.sh` script.

A failure in a test does not immediately mean that a problem is present, especially when results are compared with respect to a reference. The tests try to check if two results are reasonably close. It is in fact difficult to compare two simulated results for the electric fields. Some small variations may be normal and due to the truncation to a given floating point precision.

## User manual

cRCWA is described in an user manual available here: https://github.com/cRCWA/cRCWA/blob/main/manual/manual/afmm.pdf

## jOptiEx

cRCWA may generate large collections of files in some situations. For instance, when exporting the modes of a multimodal waveguide or structure. JOptiEx is a Java tool that allows to explore very rapidly a collection of files and represent them rapidly. It is not meant as a full-fledged scientific representation tool, but it is more a tool to explore the results to select those to be represented more carefully.

jOptiEx is a Java program. To run it, you need to have a JRE available in your computer and type `java -jar jOptiEx.jar` from the directory where `jOptiEx.jar` is present.


## Python integration

cRCWA can be used as a Python module 