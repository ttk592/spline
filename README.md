## C++ cubic spline interpolation

![cubic C2 spline](https://kluge.in-chemnitz.de/opensource/spline/cubic_c2_spline_git.png)

This is a lightweight implementation of **cubic splines**
to interpolate points f(x<sub>i</sub>) = y<sub>i</sub> with
the following features.

* available spline types:
  * **cubic C<sup>2</sup> splines**: global, twice continuously differentiable 
  * **cubic Hermite splines**: local, continuously differentiable (C<sup>1</sup>)
* boundary conditions: **first** and **second order** derivatives can be specified, **not-a-knot** condition, periodic condition is not implemented
* extrapolation
  * linear: if first order derivatives are specified or 2nd order = 0
  * quadratic: if 2nd order derivatives not equal to zero specified
* monotonicity can be enforced (when input is monotonic as well)

### Usage
The library is a header-only file with no external dependencies and can
be used like this:

```C++
#include <vector>
#include "spline.h"
...
    std::vector<double> X, Y;
    ...
    // default cubic spline (C^2) with natural boundary conditions (f''=0)
    tk::spline s(X,Y);			// X needs to be strictly increasing
    double value=s(1.3);		// interpolated value at 1.3
    double deriv=s.deriv(1,1.3);	// 1st order derivative at 1.3
    std::vector<double> solutions = s.solve(0.0);	// solves s(x)=0.0
    ...
```

The constructor can take more arguments to define the spline, e.g.:
```C++
    // cubic Hermite splines (C^1) with enforced monotonicity and
    // left curvature equal to 0.0 and right slope equal 1.0
    tk::spline s(X,Y,tk::spline::cspline_hermite, true,
                 tk::spline::second_deriv, 0.0,
                 tk::spline::first_deriv, 1.0);
```
This is identical to (must be called in that order):
```C++
    tk::spline s;
    s.set_boundary(tk::spline::second_deriv, 0.0,
                   tk::spline::first_deriv, 1.0);
    s.set_points(X,Y);
    s.make_monotonic();
```


### Spline types
Splines are piecewise polynomial functions to interpolate points
(x<sub>i</sub>, y<sub>i</sub>). In particular, cubic splines can
be represented as
* f(x) = a<sub>i</sub> + b<sub>i</sub> (x-x<sub>i</sub>) + c<sub>i</sub> (x-x<sub>i</sub>)<sup>2</sup> + d<sub>i</sub> (x-x<sub>i</sub>)<sup>3</sup>, for all x in [x<sub>i</sub>,  x<sub>i+1</sub>)
* f(x<sub>i</sub>)=y<sub>i</sub>

The following splines are available.

* `tk::spline::cspline`: cubic C<sup>2</sup> spline
  * twice continuously differentiable, e.g. f'(x<sub>i</sub>) and f''(x<sub>i</sub>) exist
  * this, together with boundary conditions uniquely determines the spline
  * requires solving a sparse equation system
  * is a global spline in the sense that changing an input point will impact the spline everywhere
  * setting first order derivatives at the boundary will break C<sup>2</sup> at the boundary
* `tk::spline::cspline_hermite`: cubic Hermite spline
  * once continuously differentiable (C<sup>1</sup>)
  * first order derivatives are specified by finite differences, e.g. on a uniform x-grid:
    * f'(x<sub>i</sub>) = (y<sub>i+1</sub>-y<sub>i-1</sub>)/(x<sub>i+1</sub>-x<sub>i-1</sub>)
  * is a local spline in the sense that changing grid points will only affect the spline around that grid point in a few adjacent segments

A function to enforce monotonicity is available as well:
* `tk::spline::make_monotonic()`: will make the spline monotonic if input grid points are monotonic
  * this function can only be called after `set_points(...)` has been called
  * it will break C<sup>2</sup> if the original spline was C<sup>2</sup> and not already monotonic
  * it will break boundary conditions if it was not monotonic in the first or last segment

### References
https://kluge.in-chemnitz.de/opensource/spline/
