#include <cstdio>
#include <cstdlib>
#include <vector>
#include "interpolation.h" // alglib

int main(int argc, char** argv)
{
    if(argc<1) {
        printf("usage: %s <>\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    std::vector<double> X(5), Y(5);
    X[0]=0.1;
    X[1]=0.4;
    X[2]=1.2;
    X[3]=1.8;
    X[4]=2.0;
    Y[0]=0.1;
    Y[1]=0.7;
    Y[2]=0.6;
    Y[3]=1.1;
    Y[4]=0.9;

    alglib::real_1d_array AX, AY;
    AX.setcontent(X.size(), &(X[0]));
    AY.setcontent(Y.size(), &(Y[0]));

    alglib::spline1dinterpolant spline;
    alglib::spline1dbuildcubic(AX, AY, X.size(), 2,0.0,2,0.0, spline);
    //alglib::spline1dbuildcubic(AX, AY, spline);


    for(size_t i=0; i<X.size(); i++) {
        printf("%f %f\n", X[i], Y[i]);
    }

    printf("\n");
    for(int i=-50; i<250; i++) {
        double x=0.01*i;
        printf("%f %f\n", x, alglib::spline1dcalc(spline,x));
    }
    printf("\n");

    return EXIT_SUCCESS;
}


