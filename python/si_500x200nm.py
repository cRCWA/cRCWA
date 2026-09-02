# The following line imports the AFMM module that should have
# been correctly installed on your machine:
import pycRCWA as crcwa

import plot2Dascii as pl


# Show the program AFMM banner and credits
crcwa.banner()
idx = crcwa.create()


#crcwa.wants(idx, "propagation")
# AFMM commands are mapped directly into Python functions:
crcwa.size(idx, 1.5e-6,1.5e-6)
crcwa.harmonics(idx, 18,18)
crcwa.wavelength(idx, 1.55e-6)


crcwa.section(idx, 2.5e-6)
# Each time in an AFMM script command there is a complex number
# to specify this is done by means of the real and imaginary part.
# In the Python access, this is handled directly by means of 
# complex variables, as in the following command:
crcwa.substrate(idx, 1.44+0j)


# Commands that are not yet accessible via Python can be accessed
# by means of 'parsescript'. You can even process a whole AFMM
# script contained in a Python string using this technique.
crcwa.parsescript(idx, "matdev la 0.0")
crcwa.pml_transf(idx, .2e-6,.2e-6,.5-0.5j)
crcwa.rectangle(idx, 3.5+0j*0,500e-9,200e-9,0,0)

# Get the refractive index distribution
struct = crcwa.inpstruct(idx, 30,25,"im")

# Here we represent the structure in the text terminal (quite crudely, but
# it gives an idea, still).
print
pl.printmap(struct)
# crcwa.bend(idx, 2e-6)
# crcwa.order(idx, 19, 20)

# Some commands give back a return value.
neff = crcwa.solve(idx)

modelistEx = crcwa.outgmodes(idx, "Ex",50,21)
modelistEy = crcwa.outgmodes(idx, "Ey",50,21)

print
k=0
for mode in modelistEx:
    print ("Mode: ",k,"n_eff=",neff[k]," |Ex| ")
    pl.printmap(modelistEx[k])
    print ("Mode: ",k,"n_eff=",neff[k]," |Ey| ")
    pl.printmap(modelistEy[k])
    k+=1
