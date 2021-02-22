// test that we can include spline.h in more than one compilation unit
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "spline.h"     // not needed here but we include it anyway
#include "myfunc.h"


int main(int argc, char** argv)
{
    double x=0.0;
    if(argc>1)
        x = atof(argv[1]);

    printf("f1(%f)=%f\n", x, myfunc1(x));
    printf("f2(%f)=%f\n", x, myfunc2(x));
    return EXIT_SUCCESS;
}

