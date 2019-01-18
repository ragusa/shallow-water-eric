# this is the gnuplot file for the seawall problem.
# note that you can execute the gnu_plot_runup.sh
# file with ./gnu_plot_runup.sh "name of EXPE txt file"
# as long as the file is in ./EXPE/
# that is, want to do something like
# 
# ./gnu_plot_seawall.sh t-20.txt
#
# the sh file will look in the correct directory
# eps
#set terminal postscript eps size 3.5,2.62 enhanced color font 'Helvetica,14' lw 3

# this sets up the general aesthetic things
set xrange [0:12]
set yrange [-0.1:0.2]
set size ratio 0.4
set key left top
set ytics -0.1,0.1,0.2
set key font ",14"
set xtics font ",14" 
set ytics font ",14" 

# this is where h + z, z get plotted with appropriate colors etc.
#sprintf("t = %0.2ss",substr(file_name,10,11)),\

#set arrow from 10.725, graph 0 to 10.725, graph 1 nohead
set style fill pattern 4 bo
plot "HplusZ.plt" every ::9 with lines lt rgb "blue" linewidth 3 title "hPz",\
     "bath.plt" every ::9 with filledcurves y1=-0.2 lt rgb "black" notitle,\
     "gauge_locations.txt" with points pt 7 ps 1 lt rgb "black" notitle

# need to set the terminal?
#set terminal postscript eps color
#set output 'GN-runup-'.sprintf("%0.2s",substr(file_name,10,11)).'.eps'
#replot

# for saving plots for GN model
#set output 'GN-runup-'.sprintf("%0.2s",substr(file_name,10,11)).'.eps'

# for saving plots for SW model
#set output 'SW-runup-'.sprintf("%0.2s",substr(file_name,10,11)).'.eps'