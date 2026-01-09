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
#   Davide Bucci, 2008-2026
#   Jérôme Michallon, 2012-2014
# 


.PHONY: src

src:
	 $(MAKE) -C src $(MAKECMDGOALS)