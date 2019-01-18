# this is the gnuplot file for the seawall problem.

# this sets up the general aesthetic things
set xrange [3:12]
set yrange [-0.2:0.15]
set size ratio 0.33
set key left top
set ytics -0.2,0.1,0.15
set key font ",16"
set xtics font ",16" 
set ytics font ",16" 

# this is where h + z, z get plotted with appropriate colors etc.
set style fill pattern 4 bo
plot "hpluszinit.plt" every ::9 with lines lt rgb "blue" linewidth 3 title "t=0s",\
     "bath.plt" every ::9 with filledcurves y1=-0.2 lt rgb "black" notitle,\
     "gauge_locations.txt" with points pt 7 ps 1 lt rgb "black" notitle


# need to set the terminal?
set terminal postscript eps color
set output 'seawall-setup.eps'
replot
