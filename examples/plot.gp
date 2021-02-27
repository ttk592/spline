# gnuplot script file

#set terminal x11 noreplotonresize
#set terminal qt size 800,800
#set terminal pngcairo size 500,500 enhanced font 'Verdana,9'
#set output "plot.png"
#set colors classic

set multiplot 
set xzeroaxis
set yzeroaxis
set key top right
set lmargin at screen 0.05	# align all plots to the left
set bmargin 0.6			# smaller bottom margins
set format x ""			# no numbers on the x-axis


set size 1.0,0.4
set origin 0.0,0.6
plot	"plot.csv" every :::0::0 title "data points" with p ps 1 pt 7 lt 1, \
	"plot.csv" every :::1::1 title "spline" with l lt 2

set size 1.0,0.2
set tmargin 0.0
set origin 0.0,0.4
plot	"plot.csv" every :::1::1 using 1:3 title "1st derivative" with l lt 2,\
	"plot.csv" every :::1::1 using 1:6 notitle with l lt 2 # should be same

set origin 0.0,0.2
plot	"plot.csv" every :::1::1 using 1:4 title "2nd derivative" with l lt 2,\
	"plot.csv" every :::1::1 using 1:7 notitle with l lt 2 # should be same


set origin 0.0,0.0
set bmargin 2.0
set format x "%g"	# numbers on the x-axis
plot	"plot.csv" every :::1::1 using 1:5 title "3rd derivative" with l lt 2

unset multiplot
pause -1
