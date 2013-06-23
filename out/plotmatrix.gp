#!/usr/env gnuplot

# plot matrix file. Invoke with:
#   $ gnuplot -e "matrixfile='file.txt'" plotmatrix.gp

set terminal pngcairo
set output matrixfile.".png"

set size ratio -1

set title "LBM D3Q19"
set xlabel "x"
set ylabel "z"
set cblabel "density"

set cbrange [0.99999:]

# Without interpolation
set pm3d map

# With interpolation
#set pm3d map interpolate 0,0


#set palette gray
splot matrixfile matrix

