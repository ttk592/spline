#!/bin/bash
./plot_avg_preserv $@ > plot.csv
return_value=$?
if [ $return_value -eq 0 ] ; then
	gnuplot plot_avg_preserv.gp
fi
