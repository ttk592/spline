#!/bin/bash

./interpol > plot.csv
./interpol_alglib >> plot.csv
gnuplot plot.txt
