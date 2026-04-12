# cRCWA

cRCWA is an efficient implementation of the Rigorous Coupled Wave Analysis (RCWA) and Aperiodic Fourier Modal Method (AFMM) techniques. The project started in 2008 with the goal of simulating a relatively large set of photonic systems and has been made open source in 2026. It is developed at the Centre de Radiofréquences, Optique et Micro-nanoélectronique des Alpes (CROMA) in Grenoble, France it has been successfully used to study diverse structures, ranging from microstructured solar cells to curved waveguides.

This open source project is distributed under the GPL v.3 license.

# Why cRCWA is different from other implementations of RCWA

The RCWA approach is a powerful method that is traditionally applied to periodic structures like Bragg gratings, photonic crystals and so on. During its development in the CROMA laboratory, cRCWA has been applied for many years to study microstructured solar cells by J. Michallon et al [1].

P. Lalanne and E. Silberstein demonstrated that RCWA can be conveniently applied to non-periodic structures [2]. They called this variant the Aperiodic Fourier Modal Method, based on the artificial periodization of the calculation window and on the introduction of Perfectly Matched Layers (PML). This allows to avoid the interaction of the structure with its replicas introduced by the periodization. 

cRCWA has been written from the beginning to implement effectively the AFMM method and to study waveguide problems, along with the more traditional periodic structures. Therefore, it should be simple to calculate the propagation modes of a waveguide of an arbitrary section, propagate the field in structures whose cross-section varies, etc. Due to the nature of the Fourier series, cRCWA is very efficient to handle cross sections of diffused waveguides for technologies such as ion-exchange in glass and lithium niobate. Furthermore, cRCWA contains implementation of state-of-the art PML with coordinate transform and normal field, to ensure the best possible performances in any situation.

cRCWA has also been the test implementation for an original development of the AFMM in cylindrical coordinates, introduced by D. Bucci, B. Martin and A. Morand in [3].

# Compile and install cRCWA

There are two ways to install cRCWA: quick installation with pre-compiled mathematical libraries or recompiling libraries on your system to improve performances. We start by presenting the quick option, but if you experience problems or if you need to optimize the performances you may explore the longer route.

## Quick install
### Quick install on Ubuntu or other Debian based Linux system

On Debian linux system, cRCWA can be quickly installed with:
~~~~
% sudo apt-get install python3 libfftw3-dev libblas-dev liblapack-dev
% ./configure
% make
% ./bin/crcwa
~~~~
cRCWA should run, if not refer to 'Troubleshooting section below'.
To test if everything is OK, you may run the automated tests:
~~~~
% cd test
%./run_all_tests.sh
~~~~
If a test fails, read the 'Self tests' section below.

### Quick install on Fodora, CentOs or other RedHat based systems

On such systems, cRCWA can be quickly installed with:
~~~~
% sudo yum install python3-devel  fftw-devel.x86_64 atlas-devel.x86_64 lapack-devel.x86_64
% ./configure
% make
% ./bin/crcwa
~~~~
cRCWA should run, if not refer to 'Troubleshooting section below'.
To test if everything is OK, you may run the automated tests:
~~~~
% cd test
%./run_all_tests.sh
~~~~
If a test fails, read the 'Self tests' section below.

### Quick install on on macOS

You can use macports to install the LAPACK, BLAS and FFTW3 libraries. The macports tool installs libraries in `/opt/local/lib`. The LAPACK library is installed in a subdirectory of this folder:
~~~~
% sudo port install lapack
% sudo port install fftw-3
% make 
~~~~
cRCWA should run, if not refer to 'Troubleshooting section below'.
To test if everything is OK, you may run the automated tests:
~~~~
% cd test
%./run_all_tests.sh
~~~~
If a test fails, read the 'Self tests' section below.
### Troubleshooting
If you were not able to build cRCWA or the python library, check that the libraries LAPACK, BLAS and FFTW3 are well installed on your system and re-lauch the configure tool. It the problem is still here, it perhabs because the libraries were not installed in usual folder. In this case, replace the config.inc file by the one corresponding to your system (config.inc.debian, config.inc.redhat, config.inc.macos). Edit the file and add manually the path to the libraries see 'detailed build instruction on Unix system below'

### Quick installation on Windows

Unfortunatly, there is no quick installation on Windows as the pre-build library of LAPACK (https://icl.utk.edu/lapack-for-windows/lapack/) lacks of gfortan library (namely libgfortran-3.dll) that we could not find. If anyone knows the work around, please let us know via GitHub's discussions :).
Instead, you need to jump to the detailed build instruction for windows (note that the prebuild of fftw3 can be directly used and requires libfftw3-3.dll to be present in the same folder of cRCWA.exe of pycRCWA.pyd). Please also read Detailed build instruction on Unix system as it gives a good information on how to optain the best computation performances.


## Detailed build instructions

cRCWA is written in C++ and can be compiled on a Unix system using GNU Make and gcc. On a practical standpoint, it has been intensively tested on Linux and on macOS. Recently Windows versions have been built, see the section 'Detailed build instruction on Window' below.

You will need the following libraries installed in your system:
- LAPACK, the classic library for linear algebra
- BLAS, the basic linear algebra package, required by LAPACK
- FFTW3, the library to efficiently compute fast Fourier transforms

These libraries can either be installed from already available pre-compiled packages (see the linux command on 'Quick installation' sections) or built from the source code. Please refer to 'Build dependency section' below for step-by-step guide for building those libraries.  

The performances of cRCWA strongly depend on the performances of the BLAS library installed in your system. Therefore, if you want to obtain the fastest calculations possible, make sure you have an optimized version of BLAS in your system. LAPACK ships with a base implementation of BLAS will work fine for test purposes, but will not deliver the best possible performance.  It is not uncommon to cut calculation times by a factor of 2 or 3 shifting from a basic BLAS implementation to an optimized one. Popular and convenient choices may be OpenBLAS or ATLAS.

### Detailed build instructions on a Unix system

Default configuration files are provided for the 3 system: debian, redhat and macOS. Copy the one corresponding to your system into confing.ind
~~~~
% cp config.inc.debian config.inc # for debian systems
% cp config.inc.redhat config.inc # for redhat systems
% cp config.inc.macos config.inc # for macOS system
~~~~
And Edit the file.

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

To create the Python module,  cRCWA also requires to know where the Python development file `Python.h` is present in your system, as well as the . This may require you to install packages such as `python...-devel` or `python...-dev`, depending on the details of your operating system.

~~~~
# This is where the Python include file Python.h is present:
PYTHONINCDIR =  [...] /3.15/include/python3.15/

# This is the Python library to be included. In macOS the file can be 
# called like libpython3.15.dylib (for a dynamically linked library).
# If you use, let's say, Python 3.11, you should change the version
# accordingly. cRCWA requires Python 3 in any case.

PYTHONLIB = -L/opt/local/Library/ [...]/3.15/lib/ -lpython3.15
~~~~

When the makefile is correctly configured, you should be able to compile cRCWA by typing `make`, while in the project main subdirectory:

~~~~
% make
~~~~

If the compile is successful, the executables are available in the `bin` directory. The Python library is put in the `python` directory.

### Detailed build instructions on a Windows

cRCWA has been successfully built on Windows using mingw and manually building LAPACK library. This section guides you though the steps.
#### Installing mingw:
Here are the installation steps:
- Download w64devkit from https://github.com/skeeto/w64devkit/releases/
- Double click on it and select where to install w64devkit
- Add the location to Windows path: 
1- Press Windows + R, type sysdm.cpl, and press Enter.
2- In the System Properties window, click on Advanced tab -> Environment Variables.
3- Under the “System variables” section, find and select Path, then click Edit.
4- Click New, then paste the path to main folder (the one provided during installation) do the same for the bin folder
5- Click OK on all open windows to save the changes.
6- Restart PC.

>Please refer to this article for more details: https://dev.to/prastha/install-mingw-w64-on-windows-11-2025-november-acg

#### Build or download pre-compiled libraries:
On windows, LAPACK library should be built as the pre-build library of LAPACK (https://icl.utk.edu/lapack-for-windows/lapack/) lacks of gfortan library (namely libgfortran-3.dll) that we could not find. Please refer to the following section 'Detailed build instructions for dependencies' to build it.
Once build, create a new 'lib' folder in cRCWA main directory and copy liblapack.a into it. If you are not planning on building OPENBLAS, also copy librefblas.a into 'lib'

Compiling FFTW3 is optional, instead you can simply download pre-compile libraries from https://www.fftw.org/install/windows.html
If you choose this option, unzip the downloaded pre-compiled library into cCRWA/lib directory and rename libfftw3-3.dll to as libfftw3.dll. Please not that with this option, you need to have libfftw3.dll present is in the same directory as cRCWA.exe or pycRCWA.pyd to execute the programs.
If you want to include them as a static library into cRCWA, refer to the following section step-by-step guide.

Compiling OPENBLAS is optional as the BLAS implementation from LAPACK could be used but computation performances will be reduced. For best computation performances, the BLAS implementation of ATLAS or OPENBLAS  should be prefered. If you skip building OPENBLAS, simply copy librefblas.a created when building LAPACK into cRCWA/lib. For details instruction on building OPENBLAS, refer to the following sections.

At this stage, the cRCWA/lib should contain:
~~~~
lib
├── liblapack.a
├── librefblas.a   
└── fftw-3.3.5-dll64
    ├── fftw3.h
    └── libfftw3.dll
~~~~
when not considering optional libraries or
~~~~
lib
├── liblapack.a
├── libopenblas.a 
└── fftw-3.3.5-dll64?
    ├── fftw3.h?
    └── libfftw3.dll?
~~~~
When using only static libraries

Then, copy config.inc.mingw into config.inc, make sure that:'LIBBLAS = -lrefblas' if you are using BLAS created by LAPACK or 'LIBBLAS = -lopenblas' is you are using BLAS created by OPENBLAS.
Then make:
~~~~
% make
~~~~
Make sure that dll libraries used during compilation are in the same folder than cRCWA. cRWCA should run!


### Detailed build instructions for dependencies on both Unix and Windows systems

#### Building lapack
Here are the download and build steps:
- Make sure that python is installed on your computer as it is used to test the library
- Download the latest version of Lapack from https://netlib.org/lapack/
- Untar.gz
- Copy make.inc.example and run make
~~~~
% cp make.inc.example make.inc
% make
~~~~
- In cRCWA main folder, create a lib folder and copy the built liblapack.a into it.

> Under windows, if you have various compilers installed, you can force to use mingw by changing CC=x86_64-w64-mingw32-gcc and FC=x86_64-w64-mingw32-gfortran
> Compilation may display an error during testing with Python when building on Windows but it should be fine to use in the next steps.

#### Building fftw3 (optional)
Here are the download and build steps:
- Download the source code from https://www.fftw.org/download.html
- Untar.gz the source code and run:

It seems that this can only be done under Linux...
~~~~
./configure --with-our-malloc16 --with-windows-f77-mangling --enable-threads --with-combined-threads --enable-portable-binary --enable-sse2 --with-incoming-stack-boundary=2 --host=x86_64-w64-mingw32
make
~~~~ 

> For more details, please refer to: https://www.fftw.org/install/windows.html

### Building OpenBlas (optional)

- Download the latest version of OpenBlas source code from https://github.com/OpenMathLib/OpenBLAS/releases
- Unzip and build:
~~~~
make
~~~~ 
> Note that building OpenBlas is very time long. Be patient...
- Copy libopenblas_haswellp-r0.3.32.a in the cRCWA/lib folder and rename it as libopenblas.a



## Self tests

cRCWA contains a certain number of self tests that can be run in the `test` directory. They can be useful to check for problems in your install. We encourage to run them often via the `run_all_tests.sh` script.

A failure in a test does not immediately mean that a problem is present, especially when results are compared with respect to a reference. The tests try to check if two results are reasonably close. It is in fact difficult to compare two simulated results for the electric fields. Some small variations may be normal and due to the truncation to a given floating point precision. If a test fails, check by hand if the result is indeed acceptable.

## How to start using cRCWA

### Simulation script

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

### Interactive program

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

### Python module

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

cRCWA is described in an user manual available here: https://github.com/cRCWA/cRCWA/blob/main/manual/LaTeX/afmm.pdf

## Tools associated to cRCWA

### jOptiEx

cRCWA may generate large collections of files in some situations. For instance, when exporting the modes of a multimodal waveguide or structure. JOptiEx is a Java tool that allows to explore very rapidly a collection of files and represent them rapidly. It is not meant as a full-fledged scientific representation tool, but it is more a tool to explore the results to select those to be represented more carefully.

jOptiEx is a Java program. To run it, you need to have a JRE available in your computer and type `java -jar jOptiEx.jar` from the directory where `jOptiEx.jar` is present.

### Slice

The `slice` utility is a tool that allows to select some lines in the output of a field propagation, contained in a file. It can be used to select a single cross-section, for instance, or reduce the size of a calculated volume:

~~~~
./slice [-lx] [-ly] [-lz] {xy|xz|yz} quota tolerance input_file output_file
~~~~

Refer to the user manual for more information.

### o2g

This utility converts Optiwave-style files (`.rid` or `.f3d`) into files that can be read with a tool like Gnuplot. Refer to the user manual for more information.

## How to contribute?

If you want to contribute to the cRCWA project, your help will be more than welcome. Here are some hints about things that can be done:
- Cite the cRCWA project in your papers.
- Benchmark cRCWA against other simulation tools.
- Improve the documentation (install instructions, use cases, ...)
- Improve the automatic test suite.
- Find and solve bugs.
- Add new functions.

If you have something in mind, do not hesitate to contact us using GitHub's discussions and describe what you plan to do and why. If you have found a bug, do not hesitate to open a new issue, by providing the minimum amount of information necessary to reproduce the bug (simulation script, etc.). Please do not provide scripts that are too long or too complicated, unless it is absolutely necessary.

If you plan to work on the code, the indentation style is K&R. We tend to avoid lines longer than 80 characters to improve readability, with the notable exception of the `commands::c_help` function.

The development is organized around three type of branches:

- `main` contains stable and tested code. pushing to main is deactivated. 
- `dev` contains ongoing development code and modifications that would soon be part of official release in `main` branch. Small modifications and bugfixes are done there. 
- `feature` branches contain added features or debugging that is under active development. The modifications here would not qualify as "small" and may require to be tested further before being merged to the `dev` branch. Once the developement is finished and tested, the branches are merged to the `dev` and deleted. 

## Bibliography

[1] - Michallon, J., Bucci, D., Morand, A., Zanuccoli, M., Consonni, V., & Kaminski-Cachopo, A. (2014). *Light trapping in ZnO nanowire arrays covered with an absorbing shell for solar cells*. Optics express, 22(S4), A1174-A1189.

[2] - Lalanne, P., & Silberstein, E. (2000). *Fourier-modal methods applied to waveguide computational problems*. Optics Letters, 25(15), 1092-1094

[3] - Bucci, D., Martin, B., & Morand, A. (2012). *Application of the three-dimensional aperiodic Fourier modal method using arc elements in curvilinear coordinates*. Optical Society of America. Journal A: Optics, Image Science, and Vision, 29(3), 367-373.
