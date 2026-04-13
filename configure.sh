# ***************************************************************************
# * automated script to configure compilation on a UNIX system              *
# *                                                                         *
# * This program detects Lapack, Blas and fftw3  			                *
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
#   Jérôme Michallon, 2012-2014 2026


#!/bin/bash
echo "Operating System Information"
echo "============================"
 
# ***************************************************************************
#                              DETECT UNIX OS
# ***************************************************************************
if [[ "$OSTYPE" == "linux-gnu"* ]]; then
 echo "Operating System: Linux"
 if [ -f /etc/os-release ]; then
  . /etc/os-release
  echo "Distribution: $PRETTY_NAME"
 fi
 echo "Kernel: $(uname -r)"
elif [[ "$OSTYPE" == "darwin"* ]]; then
 echo "Operating System: macOS"
 echo "Version: $(sw_vers -productVersion)"
echo "Build: $(sw_vers -buildVersion)"
else
  echo "Unknown OS: $OSTYPE"
fi
 
echo "Architecture: $(uname -m)"
echo "============================"

IsDebian=false
IsRedHat=false
IsMacOS=false
if [[ $PRETTY_NAME == *"Debian"* ]]; then
  echo "Debian detected!"
  IsDebian=true
elif [[ $PRETTY_NAME == *"Ubuntu"* ]]; then
  echo "Ubuntu detected!"
  IsDebian=true
elif [[ $PRETTY_NAME == *"Fedora"* ]]; then
  echo "Ubuntu detected!"
  IsRedHat=true
elif [[ $PRETTY_NAME == *"RedHat"* ]]; then
  echo "RedHat detected!"
  IsRedHat=true
elif [[ "$OSTYPE" == "darwin"* ]]; then
  echo "MacOS detected!"
  IsMacOS=true
fi

# ***************************************************************************
#                       DEFINE THE SEARCH FUNCTION :
# ***************************************************************************

# unix function to get path to libs
find_libs() {
	# pathToLib=$(dpkg -S "$@" | awk -F ": " '{print $2; exit}')
	pathToLib=$(find /usr -name "$@" 2>/dev/null | head -1)
	echo $(dirname $pathToLib)
}

# ***************************************************************************
#                         SPECIFIC OS PARAMETERS :
# ***************************************************************************


# ***************************************************************************
#                         GENERIC PARAMETERS :
# ***************************************************************************

# Tested and validated on Debian and RedHat :)
LIBLAPACKDIR=$(find_libs liblapack.so)
LIBBLASDIR=$(find_libs libblas.so)
LIBFFTW3DIR=$(find_libs libfftw3.so)
FFTW3INCDIR=$(find_libs fftw3.h)

PYVERSION=$(python3 -V | awk -F "[. ]" '{print $2 "." $3; exit}')
LIBPYTHONDIR="/usr/lib/python$PYVERSION"
PYTHONINCDIR=$(find_libs Python.h)
LIBPYTHON=-lpython$PYVERSION

# ***************************************************************************
#                         OUTPUT TO USER :
# ***************************************************************************
error=false
if [[ $LIBLAPACKDIR == "" ]]; then
  echo "Lib LAPACK not found"
  error=true
  if [ "$IsDebian" = true ]; then
    echo "Please run 'sudo apt-get install liblapack-dev'"
  elif [ "$IsRedHat" = true ]; then
    echo "Please run 'yum install lapack-devel.x86_64'"
  elif [ "$IsMacOS" = true ]; then
    echo "Please run 'sudo port install lapack'"
  fi
else
	echo "--> lib LAPACK located."
fi

if [[ $LIBBLASDIR == "" ]]; then
  echo "Lib BLAS not found"
  error=true
  if [ "$IsDebian" = true ]; then
    echo "Please run 'sudo apt-get install libblas-dev'"
  elif [ "$IsRedHat" = true ]; then
    echo "Please run 'yum install atlas-devel.x86_64 '"
  elif [ "$IsMacOS" = true ]; then
    echo "It is supposed to be installed with lapack?!"
  fi
else
	echo "--> lib BLAS located."
fi

if [[ $LIBFFTW3DIR == "" ]]; then
  echo "Lib LIBFFTW3 not found"
  error=true
  if [ "$IsDebian" = true ]; then
    echo "Please run 'sudo apt-get install libfftw3-dev'"
  elif [ "$IsRedHat" = true ]; then
    echo "Please run 'yum install fftw-devel.x86_64'"
  elif [ "$IsMacOS" = true ]; then
    echo "Please run 'sudo port install fftw-3'"
  fi
else
	echo "--> lib LIBFFTW3 located."
fi


if [ "$LIBPYTHONDIR" == "" ] || [ "$PYTHONINCDIR" == "" ]; then
  echo "Lib LIBPYTHON not found"
  error=true
  if [ "$IsDebian" = true ]; then
    echo "Please run 'sudo apt-get install python3'"
  elif [ "$IsRedHat" = true ]; then
    echo "Please run 'yum install python3-devel '"
  elif [ "$IsMacOS" = true ]; then
    echo "Please run 'sudo port install sudo port install python314'"
  fi
else
	echo "--> lib LIBPYTHON located."
fi

if [ "$error" = false ]; then
  echo "All good."
  echo "You can run 'make'"
  error=true
fi

# ***************************************************************************
#                           WRITE THE CONFIGURE.inc FILE :
# ***************************************************************************
# create configure.inc from configure.inc.auto-base:
cp config.inc.auto-base config.inc

echo "# ==================================" >>  config.inc
echo "# Automatically added configuration" >>  config.inc
echo "# ==================================" >>  config.inc
echo "LIBLAPACKDIR=$LIBLAPACKDIR" >>  config.inc
echo "LIBBLASDIR=$LIBBLASDIR" >>  config.inc
echo "LIBFFTW3DIR=$LIBFFTW3DIR" >>  config.inc
echo "FFTW3INCDIR=$FFTW3INCDIR" >>  config.inc
echo "PYVERSION=$PYVERSION" >>  config.inc
echo "LIBPYTHONDIR=$LIBPYTHONDIR" >>  config.inc
echo "PYTHONINCDIR=$PYTHONINCDIR" >>  config.inc
echo "LIBPYTHON=$LIBPYTHON" >>  config.inc


