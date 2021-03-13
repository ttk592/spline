#include <cstdio>
#include <cstdlib>
#include <vector>
#include "spline.h"

int main(int argc, char** argv)
{
    double x=0.0;
    if(argc>1)
        x = atof(argv[1]);

    std::vector<double> X(5), Y(5);
    X[0]=-0.1; X[1]=0.4; X[2]=1.2; X[3]=1.8; X[4]=2.0;
    Y[0]=0.1; Y[1]=0.7; Y[2]=0.6; Y[3]=1.1; Y[4]=0.9;


    tk::spline s1(X,Y); // defaults to C^2 cubic spline (spline::cspline)
    tk::spline s2(X,Y,tk::spline::cspline_hermite);
    tk::spline s3(X,Y,tk::spline::linear);

    printf("spline(%.3f): cubic (C^2) = %.3f, cubic hermite = %.3f, linear = %.3f\n", x, s1(x), s2(x), s3(x));

    return EXIT_SUCCESS;
}
