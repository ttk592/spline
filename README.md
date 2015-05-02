# spline
c++ cubic spline library

http://kluge.in-chemnitz.de/opensource/spline/

It generates a piecewise polynomial function of degree 3 and is
continuously differentiable everywhere. Boundary conditions are
zero-curvature at the end points. It extrapolates linearly.


## Usage
The library is a header-only file and can be used like this:

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
