# The following line imports the AFMM module that should have
# been correctly installed on your machine:
import pycRCWA as afmm

import plot2Dascii as pl


# Show the program AFMM banner and credits
afmm.banner()

xsize = 1.5e-6
ysize = 1.5e-6

# AFMM commands are mapped directly into Python functions:
afmm.size(xsize,ysize)
afmm.harmonics(31,31)
afmm.wavelength(1.55e-6)

xwgsize=0.5e-6
ywgsize=0.2e-6


npx = 1001
npy = 1001
index = [[1.44+0j] * npy for i in range(npx)]

xstart = int(npx/2-xwgsize/xsize*npx/2)
xend = xstart+int(xwgsize/xsize*npx)

ystart = int(npy/2-ywgsize/ysize*npy/2)
yend = ystart+int(ywgsize/ysize*npy)

for x in xrange(xstart, xend):
    for y in xrange(ystart, yend):
        index[y][x]=3.5+0j

afmm.indmatrix(index)

afmm.parsescript("matdev la 0.0")

# Get the refractive index distribution
struct = afmm.inpstruct(30,25,"im")

# Here we represent the structure in the text terminal (quite crudely, but
# it gives an idea, still).
print
pl.printmap(struct)

afmm.lowindex(1.44-0.1j)
afmm.highindex(3.5+0.1j)

# Some commands give back a return value.
neff = afmm.solve()

modelistEx = afmm.outgmodes("Ex",78,21)
modelistEy = afmm.outgmodes("Ey",78,21)

print
k=0
for mode in modelistEx:
    print "Mode: ",k,"n_eff=",neff[k]," |Ex| "
    pl.printmap(modelistEx[k])
    print "Mode: ",k,"n_eff=",neff[k]," |Ey| "
    pl.printmap(modelistEy[k])
    k+=1
