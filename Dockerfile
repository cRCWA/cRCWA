# syntax=docker/dockerfile:1
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    git g++ make liblapack-dev libblas-dev libfftw3-dev \
    python3 python3-dev python3-venv nano \
    && rm -rf /var/lib/apt/lists/*

# 2. Setup standard /opt/local legacy layout expected by your custom configuration
RUN mkdir -p /opt/local/lib/lapack && \
    ln -s /usr/lib/x86_64-linux-gnu/liblapack.so /opt/local/lib/lapack/liblapack.so && \
    ln -s /usr/lib/x86_64-linux-gnu/libblas.so /opt/local/lib/lapack/libblas.so && \
    ln -s /usr/lib/x86_64-linux-gnu/libfftw3.so /opt/local/lib/libfftw3.so

WORKDIR /app

# 3. Clone the cRCWA repository via HTTPS
RUN git clone https://github.com/cRCWA/cRCWA.git .

# 5. Overwrite configure.inc with placeholders, then dynamically substitute Python version
RUN cat << 'EOF' > configure.inc
# ***************************************************************************
# CONFIGURATION SECTION
# ***************************************************************************

CC = gcc
LD = g++
SHELL = /bin/sh

# This is the directory where the library files are provided. [cite: 3]
# if the LAPACK, BLAS and FFTW3 libraries are stored in the same place in [cite: 3]
# your system. [cite: 3]
# This is very often the case.
 
OPTIONS = -D LAPACK_ADD_FINAL_UNDERSCORE -D BLAS_ADD_FINAL_UNDERSCORE

LIBLAPACKDIR = /opt/local/lib/lapack
LIBLAPACK = -llapack
LIBBLASDIR = /opt/local/lib/lapack
LIBBLAS = -lblas
LIBFFTW3DIR = /opt/local/lib
LIBTHREAD = -lpthread

# Placeholders to be substituted dynamically below
PYTHONINCDIR = /usr/include/pythonPY_VER_PLACEHOLDER
PYTHONLIB = -L/usr/lib/x86_64-linux-gnu -lpythonPY_VER_PLACEHOLDER

# ***************************************************************************

# For multi-thread operation, we need:
LIBTHREAD = -lpthread

# Customize the options following the name conventions in your library.
# If in doubt, check with the nm command on the library file installed in your 
# system. 
# Variants may require a leading or final underscore. Remember that C 
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
EOF

# NEW STEP: Detect Python version and substitute PY_VER_PLACEHOLDER in configure.inc
RUN PY_VER=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')") && \
    sed -i "s/PY_VER_PLACEHOLDER/${PY_VER}/g" configure.inc

# 6. Patch src/Makefile to enforce the correct linking order and inject $(OPTIONS)
RUN sed -i '/pycrcwa:/,+3c\\
pycrcwa:\n\t$(CXX) $(OPTIONS) -shared -fPIC -I$(PYTHONINCDIR) afmm_python_module.cpp $(OFILES) $(PYTHONLIB) \\\n\t-L$(LIBLAPACKDIR) $(LIBLAPACK) -L$(LIBBLASDIR) $(LIBBLAS) -L$(LIBFFTW3DIR) -lfftw3 $(LIBTHREAD) -o pycRCWA.so\n\tmv pycRCWA.so ../python' src/makefile

# 7. Pre-compile everything so it's completely ready to use out-of-the-box
RUN make 

RUN cd src/ && make pycrcwa

CMD ["/bin/bash"]