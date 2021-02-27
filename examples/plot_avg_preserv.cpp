// plot average preserving splines:
//      this only generates output which can then be plotted using gnuplot

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
            "\t[-m] force piece-wise positivity/negativity\n"
            "\t[-x <x1> <x2> <x3> [x4] ... [xn]]\n"
            "\t[-y <y1> <y2> [y3] ... [yn-1]] \n", progname);
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
        X.resize(Y.size()+1);
        for(size_t i=0; i<X.size(); i++) {
            X[i] = (double)i;
        }
    }
    if(x_provided==true && y_provided==false) {
        fprintf(stderr, "X vector provided but not Y\n");
        exit(EXIT_FAILURE);
    }
    // check that X and Y have the correct size
    if(X.size()!=Y.size()+1 || X.size()<3) {
        fprintf(stderr, "X needs to have one element more than Y and Y needs to have at least 2 elements\n");
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

// to generate some sample data points for a histogram
double rand_unif()
{
    return (double)std::rand()/RAND_MAX;
}
double my_rand()
{
    double sum=0.0;
    for(int i=0; i<10; i++) {
        sum+=rand_unif()-0.5;
    }
    return sum;
}
int find_closest(const std::vector<double> X, double x)
{
    std::vector<double>::const_iterator it;
    it=std::upper_bound(X.begin(),X.end(),x);       // *it > x
    int idx = int(it-X.begin())-1;                    // m_x[idx] <= x
    return idx;
}
void generate_histogram(std::vector<double>& X, std::vector<double>& avg)
{
    const int num_runs = 1000;
    const int num_buckets = 20;
    const double x_min = -3.0;
    const double x_max = 3.0;
    X.resize(num_buckets+1);
    avg.resize(num_buckets);
    for(int i=0; i<=num_buckets; i++) {
        X[i] = x_min + (double)i * (x_max-x_min)/num_buckets;
    }
    for(int i=0; i<num_buckets; i++) {
        avg[i] = 0.0;
    }

    for(int i=0; i<num_runs; i++) {
        double x = my_rand();
        int idx=find_closest(X, x);
        if(0<=idx && idx<num_buckets) {
            avg[idx]++;
        }
    }
    for(int i=0; i<num_buckets; i++) {
        avg[i]/=num_buckets;
    }
}


int main(int argc, char** argv)
{
    // default parameters
    std::vector<double> X, avg;
    generate_histogram(X,avg);
    tk::spline::spline_type type = tk::spline::cspline;
    bool make_monotonic = false;

    // override defaults with command line arguments if supplied
    parse_args(argc, argv, X, avg, type, make_monotonic);

    // build Y-vector as the integral over the averages avg
    std::vector<double> Y(avg.size()+1);
    Y[0]=0.0;
    for(size_t i=1; i<Y.size(); i++) {
        Y[i] = Y[i-1] + avg[i-1]*(X[i]-X[i-1]);
    }

    // apply spline to X-Y vectors
    tk::spline s;
    s.set_boundary(tk::spline::second_deriv, 0.0,
                   tk::spline::second_deriv, 0.0);
    s.set_points(X,Y,type);     // this calculates all spline coefficients
    if(make_monotonic) {
        s.make_monotonic();     // adjusts spline coeffs to be monotonic
    }

    // the derivative of our cubic spline s then preserves the averages avg
    // this outputs the data to be used with gnuplot
    printf("# input x, input avg\n");
    for(size_t i=0; i<avg.size(); i++) {
        printf("%f %f\n", X[i], avg[i]);  // input grid points
    }
    printf("%f %f\n", X.back(), avg.back());
    printf("\n");
    int n = 1000;    // number of grid points to plot the spline
    double margin = 0.0;
    double xmin = X[0] - margin;
    double xmax = X.back() + margin;
    printf("# x, s'(x), s''(x), s'''(x)\n");
    for(int i=0; i<n; i++) {
        double x = xmin + (double)i*(xmax-xmin)/(n-1);
        printf("%f %f %f %f\n", x, s.deriv(1,x), s.deriv(2,x), s.deriv(3,x));
    }
    printf("\n");

    return EXIT_SUCCESS;
}


