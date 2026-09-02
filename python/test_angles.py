# The following line imports the AFMM module that should have
# been correctly installed on your machine:
import pycRCWA as afmm
import math

import plot2Dascii as pl


# Show the program AFMM banner and credits
afmm.banner()
idx=afmm.create()

afmm.wants(idx,"propagation")
# AFMM commands are mapped directly into Python functions:
afmm.size(idx, 1.5e-6,1.5e-6)
afmm.harmonics(idx, 2,2)
afmm.wavelength(idx, 1.55e-6)

afmm.section(idx, 0.5e-6)
# Each time in an AFMM script command there is a complex number
# to specify this is done by means of the real and imaginary part.
# In the Python access, this is handled directly by means of 
# complex variables, as in the following command:
n0=1.44
afmm.substrate(idx, n0+0j)
thetax=45*math.pi/18
thetay=30*math.pi/18


afmm.angles(idx, n0, thetax, thetay)

# Commands that are not yet accessible via Python can be accessed
# by means of 'parsescript'. You can even process a whole AFMM
# script contained in a Python string using this technique.
afmm.parsescript(idx, "matdev la 0.0")
afmm.pml_transf(idx, .2e-6,.2e-6,.5-0.5j)
afmm.rectangle(idx, 3.5+0j*0,500e-9,200e-9,0,0)

afmm.section(idx, 0.5e-6)
afmm.substrate(idx, 1.44+0j)
afmm.parsescript(idx, "matdev la 0.0")
afmm.pml_transf(idx, .2e-6,.2e-6,.5-0.5j)
afmm.rectangle(idx, 3.5+0j*0,500e-9,150e-9,0,0)


# Get the refractive index distribution
struct = afmm.inpstruct(idx, 30,25,"im")

# Some commands give back a return value.
afmm.carpet(idx)

neff = afmm.solve(idx)
afmm.parsescript(idx, "assemble")
afmm.parsescript(idx, "excitation f cy 1 0 0 0");
afmm.parsescript(idx, "propagation Ex2 m 100e-9 101 101 test_results_x.txt");
