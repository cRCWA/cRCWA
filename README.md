# cRCWA

cRCWA is an efficient implementation of the Rigorous Coupled Wave Analysis (RCWA) and Aperiodic Fourier Modal Method (AFMM) techniques. The project started in 2008 with the goal of simulating a relatively large set of photonic systems and has been made open source in 2026. It is developed at the Centre de Radiofréquences, Optique et Micro-nanoélectronique des Alpes (CROMA) in Grenoble, France it has been successfully used to study diverse structures, ranging from microstructured solar cells to curved waveguides.

This open source project is distributed under the GPL v.3 license.

## Why cRCWA is different from other implementations of RCWA

The RCWA approach is a powerful method that is traditionally applied to periodic structures like Bragg gratings, photonic crystals and so on. During its development in the CROMA laboratory, cRCWA has been applied for many years to study microstructured solar cells by J. Michallon et al [1].

P. Lalanne and E. Silberstein demonstrated that RCWA can be conveniently applied to non-periodic structures [2]. They called this variant the Aperiodic Fourier Modal Method, based on the artificial periodization of the calculation window and on the introduction of Perfectly Matched Layers (PML). This allows to avoid the interaction of the structure with its replicas introduced by the periodization. 

cRCWA has been written from the beginning to implement effectively the AFMM method and to study waveguide problems, along with the more traditional periodic structures. Therefore, it should be simple to calculate the propagation modes of a waveguide of an arbitrary section, propagate the field in structures whose cross-section varies, etc. Due to the nature of the Fourier series, cRCWA is very efficient to handle cross sections of diffused waveguides for technologies such as ion-exchange in glass and lithium niobate. Furthermore, cRCWA contains implementation of state-of-the art PML with coordinate transform and normal field, to ensure the best possible performances in any situation.

cRCWA has also been the test implementation for an original development of the AFMM in cylindrical coordinates, introduced by D. Bucci, B. Martin and A. Morand in [3].

## Compile and install on a Unix system

cRCWA is written in C++ and can be compiled on a Unix system using GNU Make and gcc. On a practical standpoint, it has been tested on Linux and on macOS.
You will need the following libraries installed in your system:
- LAPACK, the classic library for linear algebra
- BLAS, the basic linear algebra package, required by LAPACK
- FFTW3, the library to efficiently compute fast Fourier transforms

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

## Compile cRCWA on macOS

You can use macports to install the LAPACK, BLAS and FFTW3 libraries. The macports tool installs libraries in `/opt/local/lib`. The LAPACK library is installed in a subdirectory of this folder.

## Self tests

cRCWA contains a certain number of self tests that can be run in the `test` directory. They can be useful to check for problems in your install. We encourage to run them often via the `run_all_tests.sh` script.

A failure in a test does not immediately mean that a problem is present, especially when results are compared with respect to a reference. The tests try to check if two results are reasonably close. It is in fact difficult to compare two simulated results for the electric fields. Some small variations may be normal and due to the truncation to a given floating point precision. If a test fails, check by hand if the result is indeed acceptable.

## How to start using cRCWA

You have three ways to use cRCWA. The first way is to write a script file containing the commands that describe the structure and launch the simulation. For instance, `substrate` defines the refractive index of the substrate, `waveglength` the vacuum wavelength, `size` the size of the calculation window, etc. The following script calculates the guided modes for a planar waveguide. Lines starting with the symbol `#` are ignored and can contain comments:

~~~~~
# We begin by defining the size of the calculation window
# Every size and distance is specified in meters.
# Here we choose a 2 microns by 200 nanometers window

size 2e-6 200e-9

# We then have to set the number of harmonics to be used.
# Since we just have one harmonic in the y axis, we can
# we can use a pretty high number of harmonics in the x
# direction

harmonics 201 1

# We set now the substrate refractive index

substrate 1.44 0

# And we define a rectangular waveguide.
# The center of the waveguide corresponds to the center
# of the calculation window.

rectangle 3.5 0 500e-9 200e-9 0 0

# We define the wavelength

wavelength 1550e-9

# We introduce PMLs.
# Note how the size of the y absorber is zero.
# Thanks to the periodicity, we are indeed representing
# a 1D slab waveguide.

pml 5 5 50e-9 0

matdev la 1.0

# And then, we launch the calculation!!!

solve
~~~~~

In general, the most important commands write results on files. The scripting language is primitive, but contains variables, loops, conditional execution and should be flexible enough to handle simple calculations. Notice that all sizes are given in meters.

The second way to use cRCWA is using bin/crcwa as an interactive program. The program accepts the commands on a command line and executes them. It can be useful for a first approach with the software or to test interactively certain commands. Notice for instance how cRCWA complains about an incorrect spelling of a command:

~~~~
cRCWA % bin/crcwa
 ***************************************************************************
 *      Aperiodic Fourier Modal Method full vectorial 3D propagation       *
 *                            version 1.5                                  *
 *                                                                         *
 *     Build date: Jan  9 2026                                             *
 *     Source revision: N/A                                                *
 *                                                                         *
 *     Davide Bucci, CROMA March 2008 - current                            *
 *     Jérôme Michallon, CROMA May 2012 - February 2014                    *
 *     MINATEC-Grenoble INP, 3, parvis Louis Neel                          *
 *     38016, Grenoble CEDEX, France                                       *
 *                                                                         *
 *     [email redacted]                                                    *
 *                                                                         *
 ***************************************************************************
Init semaphore OK.
Reading file: stdin (interactive mode)
0 > size 10e-6 10e-6
Calculation window size set to: 1e-05 m x 1e-05 m.
1 > wavelenght 1.55e-6
Could not recognize the command wavelenght at line 2
2 > wavelength 1.55e-6
Wavelength set to: 1.55e-06 m.
3 > size 10e-6 10e-6
Calculation window size set to: 1e-05 m x 1e-05 m.
4 > quit
Stopping the program execution
~~~~

It is obviously possible to execute a script from the interactive mode using the `load` command. If you like this way to use cRCWA, you may find useful the GNU utility `rlwrap`, if it is available on your system:

~~~~
rlwrap crcwa
~~~~

The third way to use cRCWA is as a Python module. The commands that are used in the interactive mode or in the script can be used in a Python program on a module. Instead of writing on a file the calculation results, the commands return Python objects:

~~~~
# The following line imports the crcwa module that should have
# been correctly installed on your machine:
import pycRCWA as crcwa

# Show the program crcwa banner and credits
crcwa.banner()

# crcwa commands are mapped directly into Python functions:
crcwa.size(2e-6,200e-9)
crcwa.harmonics(201,1)
crcwa.wavelength(1.55e-6)

# Each time in an crcwa script command there is a complex number
# to specify this is done by means of the real and imaginary part.
# In the Python access, this is handled directly by means of 
# complex variables, as in the following command:
crcwa.substrate(1.44+0j)

# Commands that are not yet accessible via Python can be accessed
# by means of 'parsescript'. You can even process a whole crcwa
# script contained in a Python string using this technique.
crcwa.parsescript("matdev la 1.0")
crcwa.rectangle(3.5+0j*0,500e-9,200e-9,0,0)

# Some commands give back a return value.
neff = crcwa.solve()
print ("\nHere the effective indices in a Python list:\n")
print (neff)
~~~~

Not all cRCWA commands are currently mapped into Python functions, consult the list in the user manual (appendix E).

## User manual

cRCWA is described in an user manual available here: https://github.com/cRCWA/cRCWA/blob/main/manual/manual/afmm.pdf

## jOptiEx

cRCWA may generate large collections of files in some situations. For instance, when exporting the modes of a multimodal waveguide or structure. JOptiEx is a Java tool that allows to explore very rapidly a collection of files and represent them rapidly. It is not meant as a full-fledged scientific representation tool, but it is more a tool to explore the results to select those to be represented more carefully.

jOptiEx is a Java program. To run it, you need to have a JRE available in your computer and type `java -jar jOptiEx.jar` from the directory where `jOptiEx.jar` is present.

## Bibliography

[1] - Michallon, J., Bucci, D., Morand, A., Zanuccoli, M., Consonni, V., & Kaminski-Cachopo, A. (2014). /Light trapping in ZnO nanowire arrays covered with an absorbing shell for solar cells/. Optics express, 22(S4), A1174-A1189.

[2] - Lalanne, P., & Silberstein, E. (2000). /Fourier-modal methods applied to waveguide computational problems/. Optics Letters, 25(15), 1092-1094

[3] - Bucci, D., Martin, B., & Morand, A. (2012). /Application of the three-dimensional aperiodic Fourier modal method using arc elements in curvilinear coordinates/. Optical Society of America. Journal A: Optics, Image Science, and Vision, 29(3), 367-373.