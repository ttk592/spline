#include <cstdio>
#include <cstdlib>
#include <vector>
#include "spline.h"

int main(int argc, char** argv)
{
    double x=0.0;
    if(argc>1)
        x = atof(argv[1]);

    std::vector<double> X = {-0.1, 0.4, 1.2, 1.8, 2.0}; // must be increasing
    std::vector<double> Y = {0.1, 0.7, 0.6, 1.1, 0.9};

    tk::spline s1(X,Y); // defaults to C^2 cubic spline (spline::cspline)
    tk::spline s2(X,Y,tk::spline::cspline_hermite);
    tk::spline s3(X,Y,tk::spline::linear);

    printf("spline(%.3f): cubic (C^2) = %.3f, cubic hermite = %.3f, linear = %.3f\n", x, s1(x), s2(x), s3(x));

    // solve for x so that f(x) = y
    double y=0.6;
    std::vector<double> sol;    // solutions
    printf("solutions to f(x) = %f\n", y);

    sol=s1.solve(y);
    printf("C^2 spline:         ");
    for(size_t i=0; i<sol.size(); i++)
        printf("%f\t", sol[i]);
    printf("\n");

    sol=s2.solve(y);
    printf("C^1 Hermite spline: ");
    for(size_t i=0; i<sol.size(); i++)
        printf("%f\t", sol[i]);
    printf("\n");

    sol=s3.solve(y);
    printf("linear :            ");
    for(size_t i=0; i<sol.size(); i++)
        printf("%f\t", sol[i]);
    printf("\n");

    return EXIT_SUCCESS;
}
