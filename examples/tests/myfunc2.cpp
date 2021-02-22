// test that we can include spline.h in more than one compilation unit
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "spline.h"

double myfunc2(double x)
{

    static std::vector<double> X = {0.1, 0.4, 1.2, 1.8, 2.0};
    static std::vector<double> Y = {0.1, 0.7, 0.6, 1.1, 0.9};
    tk::spline s;
    s.set_boundary(tk::spline::second_deriv, 0.0,
                   tk::spline::second_deriv, 0.1);
    s.set_points(X,Y);
    return s(x);
}
