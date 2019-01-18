# this is the gnuplot file for the run up problem.
# note that you can execute the gnu_plot_runup.sh
# file with ./gnu_plot_runup.sh "name of EXPE txt file"
# as long as the file is in ./EXPE/
# that is, want to do something like
# 
# ./gnu_plot_runup.sh t-20.txt
#
# the sh file will look in the correct directory

###################################################
#### this sets up the general aesthetic things ####

set xrange [-20:12]
set yrange [-0.2:0.6]
set size ratio 0.5
set key left top
set ytics -0.2,0.2,0.6
set key font ",19"
set xtics font ",16" 
set ytics font ",16"
###################################################


# this is where h + z, z get plotted with appropriate colors etc

set style fill pattern 4 bo
plot "HplusZ.plt" every ::9 with lines lt rgb "blue" linewidth 3 title sprintf("t^{*} = %0.2ss",substr(file_name,10,11)),\
     "bath.plt" every ::9 with filledcurves y1=-0.2 lt rgb "black" notitle,\
     file_name with points pt 7 ps 1 lt rgb "red" title "EXP"

######################################################
####            for saving files                  ####

set terminal postscript eps size 4.0,2.0 enhanced color \
    font 'Helvetica,20' linewidth 2
set output 'SW-runup-'.sprintf("%0.2s",substr(file_name,10,11)).'.eps'
replot
######################################################
