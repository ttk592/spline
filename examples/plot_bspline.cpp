// plot cubic b-splines (basis splines)

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <getopt.h>
#include "spline.h"


// very simple implementation of cubic bspline basis functions
// operator()(double x) returns linear combination:
//      coeff_0 * bspline[0,4](x) + ... + coeff_n * bspline[n, n+4](x)
//  where
//      bspline[i,i+4] is the cubic basis-spline defined for all x on
//      the interval [i,i+4] and zero outside
class cubic_bspline
{
private:
    std::vector<double> m_coeff;
    // cardinal b-spline of degree 3, i.e. is piecewise cubic polynomial
    // on the interval [0,1], [1,2], [2,3] and [3,4] and zero outside
    // and is twice cont differentiable on the whole real line [-inf, inf]
    double static cubic_bspline_basis(double x)
    {
        double val = 0.0;
        if(x<0) {
            val = 0.0;
        } else if(x<1.0) {
            double z=x;
            val = z*z*z/6.0;
        } else if(x<2.0) {
            double z=x-1.0;
            val = (-3.0*z*z*z+3.0*z*z+3.0*z+1.0)/6.0;
        } else if(x<3.0) {
            double z=x-2.0;
            val = ( 3.0*z*z*z-6.0*z*z+0.0*z+4.0)/6.0;
        } else if(x<4.0) {
            double z=x-3.0;
            val = (-1.0*z*z*z+3.0*z*z-3.0*z+1.0)/6.0;
        } else {
            val = 0.0;
        }
        return val;
    }
public:
    cubic_bspline(const std::vector<double>& coeff)
    {
        m_coeff=coeff;
    }
    double operator() (double x) const
    {
        double sum=0.0;
        for(size_t i=0; i<m_coeff.size(); i++) {
            double offset = (double) i;
            sum+=m_coeff[i]*cubic_bspline_basis(x-offset);
        }
        return sum;
    }
};


void print_usage(const char* progname)
{
    printf( "usage: %s: \n"
            "\t[-h] help\n"
            "\t[-t <spline type: cspline, hermite, linear>]\n"
            "\t[-m] force piece-wise monotonicity \n"
            "\t[-b <beta1> [beta2] ... [beta n]] bspline coefficients\n"
            , progname);
}

void parse_args(int argc, char** argv,
                std::vector<double>& coeffs,
                tk::spline::spline_type& type, bool& make_monotonic)
{
    int opt;
    while( (opt = getopt(argc, argv, "ht:mb:")) != -1) {
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
        case 'b':
            coeffs.resize(0);
            optind--;
            for( ; optind < argc && *argv[optind] != '-'; optind++) {
                coeffs.push_back(atof(argv[optind]));
            }
            break;
        default:
            // we won't actually get here
            print_usage(argv[0]);
            exit(EXIT_FAILURE);
            break;
        }
    }
    // check that coeffs have the correct size
    if(coeffs.size()<1) {
        fprintf(stderr, "bspline coefficient vector (-a) need to have at least one element\n");
        exit(EXIT_FAILURE);
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
    std::vector<double> coeffs = {1.0, 0.0, 0.0, 0.0, 0.0, 1.0};
    tk::spline::spline_type type = tk::spline::cspline;
    bool make_monotonic = false;

    // override defaults with command line arguments if supplied
    parse_args(argc, argv, coeffs, type, make_monotonic);

    // setup bspline
    cubic_bspline bs(coeffs);

    // calculate corresponding X and Y coordinates from the b-splines
    std::vector<double> X(coeffs.size()+4), Y(coeffs.size()+4);
    for(size_t i=0; i<X.size(); i++) {
        X[i] = (double) i;
        Y[i] = bs(X[i]);
    }

    // setup spline
    tk::spline s(X,Y,type);
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
    printf("# x, s(x), s'(x), s''(x), bs(x), bs'(x), bs''(x)\n");
    for(int i=0; i<n; i++) {
        double x = xmin + (double)i*(xmax-xmin)/(n-1);
        printf("%f %f %f %f %f %f %f\n", x, s(x),
               s.deriv(1,x), s.deriv(2,x), bs(x),
               fderiv(bs,1,x), fderiv(bs,2,x));
    }
    printf("\n");

    return EXIT_SUCCESS;
}


