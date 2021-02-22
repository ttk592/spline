// plot splines: generates output which can be plotted using gnuplot

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <getopt.h>
#include "spline.h"

void print_usage(const char* progname)
{
    printf( "usage: %s: \n"
            "\t[-h] help\n"
            "\t[-t <spline type: cspline, hermite, linear>]\n"
            "\t[-m] force piece-wise monotonicity \n"
            "\t[-x <x1> <x2> <x3> [x4] ... [xn]]\n"
            "\t[-y <y1> <y2> <y3> [y4] ... [yn]] \n", progname);
}

void parse_args(int argc, char** argv,
                std::vector<double>& X, std::vector<double>& Y,
                tk::spline::spline_type& type, bool& make_monotonic)
{
    int opt;
    bool x_provided=false;
    bool y_provided=false;
    while( (opt = getopt(argc, argv, "ht:mx:y:")) != -1) {
        switch( opt ) {
        case 'h':
            print_usage(argv[0]);
            exit(EXIT_SUCCESS);
            break;
        case 't':
            if(std::string(optarg)=="cspline") {
                type=tk::spline::cspline;
            } else if (std::string(optarg)=="hermite") {
                type=tk::spline::cspline_hermite;
            } else if (std::string(optarg)=="linear") {
                type=tk::spline::linear;
            } else {
                print_usage(argv[0]);
                exit(EXIT_FAILURE);
            }
            break;
        case 'm':
            make_monotonic=true;
            break;
        case 'x':
            x_provided=true;
            X.resize(0);
            optind--;
            for( ; optind < argc && *argv[optind] != '-'; optind++) {
                X.push_back(atof(argv[optind]));
            }
            break;
        case 'y':
            y_provided=true;
            Y.resize(0);
            optind--;
            for( ; optind < argc && *argv[optind] != '-'; optind++) {
                Y.push_back(atof(argv[optind]));
            }
            break;
        default:
            // we won't actually get here
            print_usage(argv[0]);
            exit(EXIT_FAILURE);
            break;
        }
    }
    // create uniform x-grid if only y-grid has been provided
    if(x_provided==false && y_provided==true) {
        X.resize(Y.size());
        for(size_t i=0; i<X.size(); i++) {
            X[i] = (double)i;
        }
    }
    if(x_provided==true && y_provided==false) {
        fprintf(stderr, "X vector provided but not Y\n");
        exit(EXIT_FAILURE);
    }
    // check that X and Y have the correct size
    if(X.size()!=Y.size() || X.size()<3) {
        fprintf(stderr, "X and Y need to have the same number of elements and at least 3 points\n");
        fprintf(stderr, "provided inputs:\n");
        fprintf(stderr, "X = {");
        for(size_t i=0; i<X.size(); i++) {
            fprintf(stderr, "%.3f, ", X[i]);
        }
        fprintf(stderr, "}\n");
        fprintf(stderr, "Y = {");
        for(size_t i=0; i<Y.size(); i++) {
            fprintf(stderr, "%.3f, ", Y[i]);
        }
        fprintf(stderr, "}\n");
        exit(EXIT_FAILURE);
    }
    // check that X is monotonic
    for(int i=1; i<(int)X.size(); i++) {
        if(X[i-1]>=X[i]) {
            fprintf(stderr, "input X needs to be strictly increasing, but");
            fprintf(stderr, " X[%i]=%.3f, X[%i]=%.3f\n", i-1, X[i-1],i,X[i]);
            exit(EXIT_FAILURE);
        }
    }
}

// finite differences for 1st and 2nd order derivatives
template<class Function>
double fderiv(const Function& f, int order, double x)
{
    if(order==1) {
        double dx=2e-8*(1.0+fabs(x));
        return (f(x+dx)-f(x-dx)) / (2.0*dx);
    } else if(order==2) {
        double dx=3e-6*(1.0+fabs(x));
        return (f(x-dx)-2.0*f(x)+f(x+dx)) / (dx*dx);
    } else {
        assert(false);
        return -1.0;
    }
}



int main(int argc, char** argv)
{
    // default parameters
    std::vector<double> X = {0.1, 0.4, 1.2, 1.8, 2.0};
    std::vector<double> Y = {0.1, 0.7, 0.75, 1.1, 0.9};
    tk::spline::spline_type type = tk::spline::cspline;
    bool make_monotonic = false;

    // override defaults with command line arguments if supplied
    parse_args(argc, argv, X, Y, type, make_monotonic);

    // setup spline
    tk::spline s;
    // boundary conditions: natural condition is 2nd deriv = 0
    // extrapolation:
    //  - at most quadratic and linear if 2nd deriv = 0
    //  - linear if 1st derivative at boundary is specified (this breaks C^2)
    s.set_boundary(tk::spline::second_deriv, 0.0,
                   tk::spline::second_deriv, 0.0);
    s.set_points(X,Y,type);     // this calculates all spline coefficients
    if(make_monotonic) {
        s.make_monotonic();     // adjusts spline coeffs to be monotonic
    }

    // evaluates spline and outputs data to be used with gnuplot
    printf("# input x, input y\n");
    for(size_t i=0; i<X.size(); i++) {
        printf("%f %f\n", X[i], Y[i]);  // input grid points
    }
    printf("\n");
    int n = 1000;    // number of grid points to plot the spline
    double xmin = X[0] - 0.5;
    double xmax = X.back() + 0.5;
    printf("# x, s(x), s'(x), s''(x), s'''(x), fdiff(s,1,x), fdiff(s,2,x)\n");
    for(int i=0; i<n; i++) {
        double x = xmin + (double)i*(xmax-xmin)/(n-1);
        printf("%f %f %f %f %f %f %f\n", x, s(x),
               s.deriv(1,x), s.deriv(2,x), s.deriv(3,x),
               fderiv(s,1,x), fderiv(s,2,x));
    }
    printf("\n");

    return EXIT_SUCCESS;
}


