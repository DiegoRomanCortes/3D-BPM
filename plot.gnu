# plots and saves to png the BPM data
#set cbrange [0:1]
set palette defined ( 0 "black", 1 "blue", 2 "cyan", 3 "green", 4 "yellow", 5 "orange", 6 "red", 7 "magenta", 8 "white" ) 

Lx = 50E-6; # width of the grid
Ly = 50E-6; # height of the grid

# call gnuplot to plot the refractive index
set pm3d map
set size ratio -1
set xrange [-Lx/2:Lx/2]
set yrange [-Ly/2:Ly/2]
set xlabel 'X (m)'
set ylabel 'Y (m)'
set title 'Refractive Index'
splot 'refractive2d.txt' using 1:2:3 with image
set term png
set output 'refractive.png'
replot
set term x11

# call gnuplot to plot the input data
set pm3d map
set xrange [-Lx/2:Lx/2]
set yrange [-Ly/2:Ly/2]
set xlabel 'X (m)'
set ylabel 'Y (m)'
set title 'Initial Gaussian Field'
splot '0.txt' using 1:2:3 with image
set term png
set output '0.png'
replot
set term x11

# call gnuplot to plot the output data
do for [num = 61:62]{
    set pm3d map
    set size ratio -1
    set xrange [-Lx/2:Lx/2]
    set yrange [-Ly/2:Ly/2]
    set xlabel 'X (m)'
    set ylabel 'Y (m)'
    set title 'Propagated Gaussian Field'
    splot ''.num .'.txt' using 1:2:3 with image
    set term png
    set output ''.num .'.png'
    replot
    set term x11
}