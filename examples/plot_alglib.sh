#!/bin/bash
./plot_alglib $@ > plot.csv
return_value=$?
if [ $return_value -eq 0 ] ; then
	gnuplot plot_alglib.gp
fi
