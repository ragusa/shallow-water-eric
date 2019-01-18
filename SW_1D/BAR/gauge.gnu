set terminal qt
set xrange [33:39]
set yrange [-0.02:0.1]
set size ratio 0.33
set key font ",19"
set xtics font ",16"
set ytics font ",16"
plot "WG1.txt" using 1:3 with lines lt rgb "blue" linewidth 3 title "h + z",\
     "./EXPE/wg1.txt" with points pt 7 ps 1 lt rgb "red" title "EXP"

set terminal qt 1
set xrange [0:10]
set yrange [-0.02:0.1]
set size ratio 0.33
set key font ",19"
set xtics font ",16"
set ytics font ",16"
plot "WG3.txt" using 1:3 with lines lt rgb "blue" linewidth 3 title "h + z",\
     "./EXPE/wg3.txt" with points pt 7 ps 1 lt rgb "red" title "EXP"

set terminal qt 2
set xrange [2:10]
set yrange [-0.02:0.1]
set size ratio 0.33
set key font ",19"
set xtics font ",16"
set ytics font ",16"
plot "WG10.txt" using 1:3 with lines lt rgb "blue" linewidth 3 title "h + z",\
     "./EXPE/wg10.txt" with points pt 7 ps 1 lt rgb "red" title "EXP"

set terminal qt 3
set xrange [0:10]
set yrange [-0.02:0.1]
set size ratio 0.33
set key font ",19"
set xtics font ",16"
set ytics font ",16"
plot "WG22.txt" using 1:3 with lines lt rgb "blue" linewidth 3 title "h + z",\
     "./EXPE/wg22.txt" with points pt 7 ps 1 lt rgb "red" title "EXP"

set terminal qt 4
set xrange [0:10]
set yrange [-0.02:0.1]
set size ratio 0.33
set key font ",19"
set xtics font ",16"
set ytics font ",16"
plot "WG28.txt" using 1:2 with lines lt rgb "blue" linewidth 3 title "h",\
     "./EXPE/wg28.txt" with points pt 7 ps 1 lt rgb "red" title "EXP"

set terminal qt 5
set xrange [0:10]
set yrange [-0.02:0.1]
set size ratio 0.33
set key font ",19"
set xtics font ",16"
set ytics font ",16"
plot "WG37.txt" using 1:2 with lines lt rgb "blue" linewidth 3 title "h",\
     "./EXPE/wg37.txt" with points pt 7 ps 1 lt rgb "red" title "EXP"

set terminal qt 4
set title "WG40"
set xrange [0:10]
set yrange [-0.02:0.1]
set size ratio 0.33
set key font ",19"
set xtics font ",16"
set ytics font ",16"
plot "WG40.txt" using 1:2 with lines lt rgb "blue" linewidth 3 title "h",\
     "./EXPE/wg40.txt" with points pt 7 ps 1 lt rgb "red" title "EXP"
