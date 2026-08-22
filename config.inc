
# ***************************************************************************
#                    UNIX BASE CONFIG FOR AUTOTOOL  
# ***************************************************************************

# Define the compilers:
CC = gcc
LD = g++
CXX = g++
SHELL = /bin/sh


LIBLAPACK = -llapack
LIBBLAS = -lblas
LIBFFTW3 = -lfftw3

# Define specific OS parameters:
EXE_EXT :=
PY_EXT := .so

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

# ***************************************************************************
# ==================================
# Automatically added configuration
# ==================================
LIBLAPACKDIR=/usr/lib/x86_64-linux-gnu/lapack
LIBBLASDIR=/usr/lib/x86_64-linux-gnu
LIBFFTW3DIR=/usr/lib/x86_64-linux-gnu
FFTW3INCDIR=/usr/include
PYVERSION=3.12
LIBPYTHONDIR=/usr/lib/python3.12
PYTHONINCDIR=/usr/include/python3.12
LIBPYTHON=-lpython3.12
