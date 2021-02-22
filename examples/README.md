## Spline example use cases

### Prerequisites
* C++11 compatible compiler
* For a few examples [ALGLIB](https://www.alglib.net/) is required

### Building
Adjust the variables `CC`, `CFLAGS` and `ALGLIB` in the `Makefile` and
execute as follows:
```
$ make 		# builds all self-contained examples
$ make more	# builds all examples which require ALGLIB
$ make tests	# builds and runs several tests
```

### Example files (self contained)

* `simple_demo.cpp`: sets up a spline and evaluates it at a single point
* `plot*.cpp`: there are several examples to plot splines in different situations.
	They all produce numerical output which can be used with [Gnuplot](http://www.gnuplot.info/) for visualisation. Shell scripts are provided to do this in one step: `./plot*.sh` (or manually: `./plot > plot.csv`, `gnuplot plot.gp`).
  * `plot.cpp`: plot splines, e.g.
	```
	$ ./plot.sh -t hermite -x 0.0 0.5 1.5 2.0 -y 1.1 0.9 0.5 0.7
	```
  * `plot_avg_preserv.cpp`: plots **average preserving** (area preserving/mass preserving) splines. These are splines which do not interpolate points but which have an average value in the interval [x<sub>i</sub>, x<sub>i+1</sub>] equal to y<sub>i</sub>.
	```
	$ ./plot_avg_preserv.sh -t cspline 
	```
  * `plot_parametric.cpp`: plots a **parametric spline**, i.e. connects points in a 2D-plane, where the x-coordinates do no longer need to be strictly increasing.
	```
	$ ./plot_parametric.sh -t hermite
	```

### Example files (which require ALGLIB)
* `bench.cpp`: performs a benchmark of creating a spline and interpolating values. E.g. on a 2.4GHz CPU, spline size 50 and 1 million spline evaluations call:
	```
	$ ./bench 2400 50 1000
	                                tk                      alglib
	random access:   loops=1e+06,   0.064s ( 153 cycl)      0.078s ( 188 cycl)
	spline creation: loops=2e+04,   0.200s (2.4e+04 cycl)   0.202s (2.4e+04 cycl)
	grid transform:  loops=2e+04,   0.250s (3.0e+04 cycl)   0.193s (2.3e+04 cycl)
	```
* `plot_alglib.cpp`: plots ALGLIB spline against this spline implementation (does not have any command line options)
	```
	$ ./plot_alglib.sh
	```

