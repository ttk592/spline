# gnuplot script file

set size ratio -1
set xzeroaxis
set yzeroaxis
set key top right
plot    "plot.csv" every :::0::0 using 1:2 notitle with p ps 2, \
        "plot.csv" every :::1::1 using 2:3 title "parametric spline" with l lt 1

pause -1
