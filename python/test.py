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
