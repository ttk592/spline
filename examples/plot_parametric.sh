#!/bin/bash
./plot_parametric $@ > plot.csv
return_value=$?
if [ $return_value -eq 0 ] ; then
	gnuplot plot_parametric.gp
fi
