
# ***************************************************************************
# CONFIGURATION SECTION
# ***************************************************************************

CC = gcc
LD = g++
SHELL = /bin/sh

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
