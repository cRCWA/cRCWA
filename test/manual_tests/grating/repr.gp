set pm3d
set view map
set terminal png size 600,400
set out "field.png"
!../../slice -lz xz 0 1e-6 propagation.txt propagation_c.txt
splot "propagation_c.txt" u 1:3:4 w pm3d title ""

