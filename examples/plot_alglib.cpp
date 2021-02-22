#include <cstdio>
#include <cstdlib>
#include <vector>
#include "interpolation.h"  // alglib
#include "spline.h"         // tk

// wrap alglib spline into a class
class alglib_spline
{
private:
    alglib::spline1dinterpolant m_spline;
public:
    alglib_spline(const std::vector<double>& X, const std::vector<double>& Y)
    {
        alglib::real_1d_array AX, AY;
        AX.setcontent(X.size(), &(X[0]));
        AY.setcontent(Y.size(), &(Y[0]));
        //alglib::spline1dbuildcubic(AX, AY, m_spline);
        alglib::spline1dbuildcubic(AX, AY, X.size(), 2,0.0,2,0.0, m_spline);
    }
    double operator()(double x) const
    {
        return alglib::spline1dcalc(m_spline,x);
    }
};

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
    if(argc<1) {
        printf("usage: %s <>\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    std::vector<double> X = {0.1, 0.4, 1.2, 1.8, 2.0};
    std::vector<double> Y = {0.1, 0.7, 0.6, 1.1, 0.9};

    alglib_spline   aspline(X,Y);
    tk::spline      tspline(X,Y);

    printf("# input x, input y\n");
    for(size_t i=0; i<X.size(); i++) {
        printf("%f %f\n", X[i], Y[i]);
    }
    printf("\n");
    printf("# x, alglib spline(x), tk spline(x), as'(x), ts'(x), as''(x), ts''(x)\n");
    for(int i=-50; i<250; i++) {
        double x=0.01*i;
        printf("%f %f %f %f %f %f %f\n", x, aspline(x), tspline(x),
               fderiv(aspline,1,x), fderiv(tspline,1,x),
               fderiv(aspline,2,x), fderiv(tspline,2,x));
    }
    printf("\n");

    return EXIT_SUCCESS;
}


