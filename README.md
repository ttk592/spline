## spline
c++ cubic spline library

http://kluge.in-chemnitz.de/opensource/spline/

It generates a piecewise polynomial function of degree 3 and is
twice continuously differentiable everywhere. Boundary conditions
default to zero-curvature at the end points. It extrapolates linearly,
if default boundary conditions are used, or otherwise extrapolation
is a quadratic function.


### Usage
The library is a header-only file and can be used like this:

```C++
#include "spline.h"
...
int main(int argc, char** argv)
{
    std::vector<double> X, Y;
    ...
    tk::spline s;
    s.set_points(X,Y);    // X needs to be sorted, strictly increasing
    double value=s(1.5);  // interpolated value at 1.5
    ...
}
```

Optionally, boundary condition can be modified via `set_boundary()`,
which needs to be called before `set_points()`.

```C++
s.set_boundary(tk::spline::second_deriv,0.0,tk::spline::first_deriv,-2.0,false);
```
