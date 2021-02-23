#!/bin/bash
./plot_bspline $@ > plot.csv
return_value=$?
if [ $return_value -eq 0 ] ; then
	gnuplot plot_bspline.gp
fi
