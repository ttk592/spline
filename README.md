## C++ cubic spline interpolation

This is a lightweight implementation of **cubic splines**
with the following features:

* available spline types:
  * cubic (classical) splines: global, twice continuously differentiable (C<sup>2</sup>)
  * cubic Hermite splines: local, continuously differentiable (C<sup>1</sup>)
* boundary conditions: first and second order derivatives can be specified
* extrapolation
  * is at most a quadratic function
  * linear when second order derivatives are set to zero at boundary points or when first derivatives are specified
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
    ...
}
```

The constructor can take more arguments to define the spline, e.g.:
```C++
    // cubic Hermite splines (C^1) with enforced monotonicity and
    // left curvature equal to 0.0 and right slope equal 1.0
    tk::spline s(X,Y,tk::spline::cspline_hermite, true,
                 tk::spline::second_deriv, 0.0,
                 tk::spline::first_deriv, 1.0);
}
```
This is identical to (must be called in that order):
```C++
    tk::spline s;
    s.set_boundary(tk::spline::second_deriv, 0.0,
                   tk::spline::first_deriv, 1.0);
    s.set_points(X,Y);
    s.make_monotonic();
}
```

### Spline types
* `tk::spline::cspline`: classical cubic spline
  * is uniquely determined by the requirement to be twice continuously differentiable (C<sup>2</sup>)
  * requires solving a sparse equation system
  * is a global spline in the sense that changing an input point will impact the spline everywhere
  * setting first order derivatives at the boundary will break C<sup>2</sup> at the boundary
* `tk::spline::cspline_hermite`: cubic hermite spline
  * is continuously differentiable (C<sup>1</sup>) and derivatives are specified on every grid point (equal to the finite differences of the 3 adjacent grid points)
  * is a local spline in the sense that changing grid points will only affect the spline around that grid point in a few adjacent segments

A function to enforce monotonicity is available as well:
* `tk::spline::make_monotonic()`: will make the spline monotonic if input grid points are monotonic
  * this function can only be called after `set_points(...)` has been called
  * it will break C<sup>2</sup> if the original spline was C<sup>2</sup> and not already monotonic
  * it will break boundary conditions if it was not monotonic in the first or last segment

### References
http://kluge.in-chemnitz.de/opensource/spline/
