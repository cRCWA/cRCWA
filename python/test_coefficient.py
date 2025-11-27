# The following line imports the AFMM module that should have
# been correctly installed on your machine:
import pyAFMM as afmm

import plot2Dascii as pl


# Show the program AFMM banner and credits
afmm.banner()


afmm.wants("propagation")
# AFMM commands are mapped directly into Python functions:
afmm.size(1.5e-6,1.5e-6)
afmm.harmonics(17,17)
afmm.wavelength(1.55e-6)

afmm.section(0.5e-6)
# Each time in an AFMM script command there is a complex number
# to specify this is done by means of the real and imaginary part.
# In the Python access, this is handled directly by means of 
# complex variables, as in the following command:
afmm.substrate(1.44+0j)


# Commands that are not yet accessible via Python can be accessed
# by means of 'parsescript'. You can even process a whole AFMM
# script contained in a Python string using this technique.
afmm.parsescript("matdev la 0.0")
afmm.pml_transf(.2e-6,.2e-6,.5-0.5j)
afmm.rectangle(3.5+0j*0,500e-9,200e-9,0,0)

afmm.section(0.5e-6)
afmm.substrate(1.44+0j)
afmm.parsescript("matdev la 0.0")
afmm.pml_transf(.2e-6,.2e-6,.5-0.5j)
afmm.rectangle(3.5+0j*0,500e-9,150e-9,0,0)


# Get the refractive index distribution
struct = afmm.inpstruct(30,25,"im")

# Some commands give back a return value.
afmm.carpet()

neff = afmm.solve()
afmm.parsescript("assemble")
afmm.parsescript("excitation f cy 1 0 0 0");
afmm.parsescript("propagation Ex2 m 100e-9 101 101 test_results_x.txt");

print(afmm.coefficient("f",1.49033));
print(afmm.coefficient_id("f",10));

