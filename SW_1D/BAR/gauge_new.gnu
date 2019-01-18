set terminal qt 3
set xrange [-2:2]
set yrange [-0.02:0.1]
set size ratio 0.33
plot "./new_update/WG1.txt" using 1:3 with lines 
replot "./EXPE/wg1.txt"

set terminal qt 4
set xrange [0:10]
set yrange [-0.02:0.1]
set size ratio 0.33
plot "./new_update/WG3.txt" using 1:3 with lines 
replot "./EXPE/wg3.txt"

set terminal qt 5
set xrange [2:10]
set yrange [-0.02:0.1]
set size ratio 0.33
plot "./new_update/WG10.txt" using 1:3 with lines
replot "./EXPE/wg10.txt"