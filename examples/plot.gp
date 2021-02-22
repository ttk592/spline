# gnuplot script file

#set terminal x11 noreplotonresize
#set terminal qt size 800,800

set multiplot 
set xzeroaxis
set yzeroaxis
set key top right
set lmargin at screen 0.05	# align all plots to the left


set size 1.0,0.4
set origin 0.0,0.6
plot	"plot.csv" every :::0::0 notitle with p ps 2, \
	"plot.csv" every :::1::1 title "spline" with l lt 1

set size 1.0,0.2
set origin 0.0,0.4
plot	"plot.csv" every :::1::1 using 1:3 title "1st derivative" with l lt 1,\
	"plot.csv" every :::1::1 using 1:6 notitle with l lt 1 # should be same

set origin 0.0,0.2
plot	"plot.csv" every :::1::1 using 1:4 title "2nd derivative" with l lt 1,\
	"plot.csv" every :::1::1 using 1:7 notitle with l lt 1 # should be same

set origin 0.0,0.0
plot	"plot.csv" every :::1::1 using 1:5 title "3rd derivative" with l lt 1

unset multiplot
pause -1
