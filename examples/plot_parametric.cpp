// plot parametric splines
//  this generates output which can be plotted using gnuplot

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
    if(x_provided!=y_provided) {
        fprintf(stderr, "need to provide both X and Y vectors or none\n");
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
    // check that two succesive points are not identical
    for(int i=1; i<(int)X.size(); i++) {
        if( X[i-1]==X[i] && Y[i-1]==Y[i]) {
            fprintf(stderr, "two successive points must not be the same, but");
            fprintf(stderr, " X[%i]=%.3f, X[%i]=%.3f, Y[%i]=%.3f, Y[%i]=%.3f\n",
                    i-1,X[i-1], i,X[i], i-1,Y[i-1], i,Y[i]);
            exit(EXIT_FAILURE);
        }
    }
}

double sqr(double x)
{
    return x*x;
}

int main(int argc, char** argv)
{
    // default parameters
    std::vector<double> X = {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 2.0};
    std::vector<double> Y = {1.0, 0.0, 0.0, 1.0, 1.0, 2.0, 2.0};
    tk::spline::spline_type type = tk::spline::cspline;
    bool make_monotonic = false;

    // override defaults with command line arguments if supplied
    parse_args(argc, argv, X, Y, type, make_monotonic);

    // setup parametric spline
    // introduce a "time variable"
    std::vector<double> T(X.size());
    T[0]=0.0;
    for(size_t i=1; i<T.size(); i++) {
        T[i] = T[i-1] + sqrt( sqr(X[i]-X[i-1]) + sqr(Y[i]-Y[i-1]) );
    }
    // define a spline for each coordinate x, y
    tk::spline sx, sy;
    sx.set_points(T,X,type);
    sy.set_points(T,Y,type);
    if(make_monotonic) {
        // adjusts spline coeffs to be piecewise monotonic where possible
        sx.make_monotonic();
        sy.make_monotonic();
    }

    // evaluates spline and outputs data to be used with gnuplot
    printf("# input x, input y\n");
    for(size_t i=0; i<X.size(); i++) {
        printf("%f %f\n", X[i], Y[i]);  // input grid points
    }
    printf("\n");
    int n = 1000;    // number of grid points to plot the spline
    double tmin = T[0] - 0.5;
    double tmax = T.back() + 0.5;
    printf("# t, sx(t), sy(t)\n");
    for(int i=0; i<n; i++) {
        double t = tmin + (double)i*(tmax-tmin)/(n-1);
        printf("%f %f %f\n", t, sx(t), sy(t));
    }
    printf("\n");

    return EXIT_SUCCESS;
}


