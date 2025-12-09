# The following line imports the AFMM module that should have
# been correctly installed on your machine:
import pycRCWA as afmm

import plot2Dascii as pl


# Show the program AFMM banner and credits
afmm.banner()


afmm.wants("propagation")
# AFMM commands are mapped directly into Python functions:
afmm.size(1.5e-6,1.5e-6)
afmm.harmonics(17,17)
afmm.wavelength(1.55e-6)


afmm.section(2.5e-6)
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

# Get the refractive index distribution
struct = afmm.inpstruct(30,25,"im")

# Here we represent the structure in the text terminal (quite crudely, but
# it gives an idea, still).
print
pl.printmap(struct)
# afmm.bend(2e-6)
# afmm.order(19, 20)

# Some commands give back a return value.
neff = afmm.solve()

modelistEx = afmm.outgmodes("Ex",50,21)
modelistEy = afmm.outgmodes("Ey",50,21)

print
k=0
for mode in modelistEx:
    print ("Mode: ",k,"n_eff=",neff[k]," |Ex| ")
    pl.printmap(modelistEx[k])
    print ("Mode: ",k,"n_eff=",neff[k]," |Ey| ")
    pl.printmap(modelistEy[k])
    k+=1
