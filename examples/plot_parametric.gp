# gnuplot script file

#set terminal pngcairo size 400,400 enhanced font 'Verdana,9'
#set output "plot.png"
#set colors classic


set size ratio -1
set xzeroaxis
set yzeroaxis
set key top right
plot    "plot.csv" every :::0::0 using 1:2 title "data points" with p ps 1 pt 7 lt 1, \
        "plot.csv" every :::1::1 using 2:3 title "parametric spline" with l lt 2

pause -1
