# gnuplot script file

#set terminal x11 noreplotonresize
#set terminal qt size 800,800

set multiplot 
set xzeroaxis
set yzeroaxis
set key top right
set lmargin at screen 0.05	# align all plots to the left


set size 1.0,0.5
set origin 0.0,0.5
plot	"plot.csv" every :::0::0 notitle with p ps 2, \
        "plot.csv" every :::1::1 using 1:5 title "b-spline" with l lt 1,\
        "plot.csv" every :::1::1 using 1:2 title "tk spline" with l lt 2


set size 1.0,0.25
set origin 0.0,0.25
plot	"plot.csv" every :::1::1 using 1:6 title "1st deriv b" with l lt 1,\
	"plot.csv" every :::1::1 using 1:3 title "1st deriv tk"  with l lt 2

set origin 0.0,0.0
plot	"plot.csv" every :::1::1 using 1:7 title "2nd deriv b" with l lt 1,\
	"plot.csv" every :::1::1 using 1:4 title "2nd deriv tk" with l lt 2

unset multiplot
pause -1
