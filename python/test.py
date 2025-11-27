# The following line imports the AFMM module that should have
# been correctly installed on your machine:
import pyAFMM as afmm

# Show the program AFMM banner and credits
afmm.banner()

# AFMM commands are mapped directly into Python functions:
afmm.size(2e-6,200e-9)
afmm.harmonics(201,1)
afmm.wavelength(1.55e-6)

# Each time in an AFMM script command there is a complex number
# to specify this is done by means of the real and imaginary part.
# In the Python access, this is handled directly by means of 
# complex variables, as in the following command:
afmm.substrate(1.44+0j)

# Commands that are not yet accessible via Python can be accessed
# by means of 'parsescript'. You can even process a whole AFMM
# script contained in a Python string using this technique.
afmm.parsescript("matdev la 1.0")
afmm.rectangle(3.5+0j*0,500e-9,200e-9,0,0)

# Some commands give back a return value.
neff = afmm.solve()
print ("\nHere the effective indices in a Python list:\n")
print (neff)
