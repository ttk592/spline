# gnuplot script file

#set term x11 noreplotonresize
#set terminal qt size 800,800


set lmargin at screen 0.05      # align all plots to the left

set multiplot

set xzeroaxis
set yzeroaxis
set grid y
set size 1.0,0.5
set origin 0.0,0.5
set key top right
set title "average preserving spline interpolation"
plot    "plot.csv" every :::0::0 using 1:2 notitle with steps lt 2, \
        "plot.csv" every :::1::1 using 1:2 title "spline" with l lt 1


set size 1.0,0.25
set origin 0.0,0.25
set notitle
plot    "plot.csv" every :::1::1 using 1:3 title "1st derivative" with l lt 1
        
set size 1.0,0.25
set origin 0.0,0.0
plot    "plot.csv" every :::1::1 using 1:4 title "2nd derivative" with l lt 1

unset multiplot

pause -1
