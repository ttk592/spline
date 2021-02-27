# gnuplot script file

#set terminal x11 noreplotonresize
#set terminal qt size 800,800
#set terminal pngcairo size 400,400 enhanced font 'Verdana,9'
#set colors classic
#set output "plot.png"

set xzeroaxis
set yzeroaxis
set key top right
set lmargin at screen 0.07	# align all plots to the left
set bmargin 0.6
set format x ""

set multiplot 

set size 1.0,0.5
set origin 0.0,0.5
set key bottom right
plot	"plot.csv" every :::0::0 title "data points" with p ps 1 lt 1 pt 7, \
        "plot.csv" every :::1::1 using 1:2 title "alglib spline" with l lt 1,\
        "plot.csv" every :::1::1 using 1:3 title "tk spline" with l lt 2


set tmargin 0
set size 1.0,0.25
set origin 0.0,0.25
set key top right
plot	"plot.csv" every :::1::1 using 1:4 title "1st deriv alglib" with l lt 1,\
	"plot.csv" every :::1::1 using 1:5 title "1st deriv tk"  with l lt 2


set format x "%.1f"
set bmargin 2.0
set origin 0.0,0.0
plot	"plot.csv" every :::1::1 using 1:6 title "2nd deriv alglib" with l lt 1,\
	"plot.csv" every :::1::1 using 1:7 title "2nd deriv tk" with l lt 2

unset multiplot
pause -1
