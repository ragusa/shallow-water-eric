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
do for [ii=1:199] {
set xrange [0:12]
set yrange [-0.1:0.2]
set size ratio 0.4
set key left top
set ytics -0.1,0.1,0.2
set key font ",14"
set xtics font ",14" 
set ytics font ",14" 

file = sprintf('hpz_movie_%0.0f.plt',ii)
set style fill solid 
plot file every ::9 with filledcurves y=-0.2 lt rgb "blue" linewidth 3 notitle, \
     "bath.plt" every ::9 with filledcurves y1=-0.2 lt rgb "black" notitle, \
     "../gauge_locations.txt" with points pt 7 ps 1 lt rgb "red" notitle
    
######################################################
####            for saving files                  ####
#eps size 4.0,2.0 enhanced color 
#set terminal png 
set term pngcairo size 1024,768
outfile = sprintf('movie.%03.0f.png',ii)
set output outfile
replot
######################################################
}
