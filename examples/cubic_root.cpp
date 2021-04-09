#include <cstdio>
#include <cstdlib>
#include <complex>
#include <vector>
#include <array>
#include <time.h>
#include <sys/time.h>
#include "spline.h"

// gsl functions
// f(x) = x^3 + a x^2 + b x + c
int gsl_poly_solve_cubic (double a, double b, double c,
                          double *x0, double *x1, double *x2);
// f(x) = a x^2 + b x + c
int gsl_poly_solve_quadratic (double a, double b, double c,
                              double *x0, double *x1);


// f(x) = a + b x + c x^2 + d x^3
std::vector<double> gsl_poly_wrapper(double a, double b, double c, double d)
{
    double x0=0.0, x1=0.0, x2=0.0;
    int num_sol=0;
    if(d==1.0) {
        num_sol = gsl_poly_solve_cubic(c,b,a,&x0,&x1,&x2);
    } else if(d!=0.0) {
        num_sol = gsl_poly_solve_cubic(c/d,b/d,a/d,&x0,&x1,&x2);
    } else { // d==0.0
        num_sol = gsl_poly_solve_quadratic(c,b,a,&x0,&x1);
    }
    std::vector<double> x;
    if(num_sol>=1)
        x.push_back(x0);
    if(num_sol>=2)
        x.push_back(x1);
    if(num_sol>=3)
        x.push_back(x2);
    return x;
}


// time of today in seconds (with micro seconds resolution)
double stoptime(void)
{
    struct timeval t;
    gettimeofday(&t,NULL);
    return (double) t.tv_sec + t.tv_usec/1000000.0;
}


// uniform distribution
double rand_unif()
{
    return (double)std::rand()/RAND_MAX;
}

// double exponential distribution
double rand_laplace(double lambda=1.0, double mu=0)
{
    double u = rand_unif();
    if(u<=0.5) {
        return mu-log(2.0*u)/lambda;
    } else {
        return mu+log(2.0*(u-0.5))/lambda;
    }
}

// round to n significant decimals
// inefficient (pow, log very slow)
double round_dec(double x, int n)
{
    double expon = floor(log10(fabs(x)));
    double scale = pow(10.0, expon-n);
    double y = round(x/scale);
    return y*scale;
}

// polynomial: f(x) = a + b*x + c*x^2 + x^3
double f(double a, double b, double c, double d, double x)
{
    return ((d*x+c)*x+b)*x+a;
}
// polynomial: f(x) = a + b*x + c*x^2 + x^3
double f(double a, double b, double c, double x)
{
    return f(a,b,c,1.0,x);
}

// calculate polynomial coefficients from exact roots
// f(x) = a + b*x + c*x^2 + x^3
void coeffs_from_roots(double& a, double&b, double&c,
                       const std::array<std::complex<double>, 3>& root)
{
    // (x-x0)*(x-x1)*(x-x2) = -x0*x1*x2 + (x0*x1+x0*x2+x1*x2)*x
    //                        -(x0+x1+x2)*x^2 + x^3
    a = -std::real(root[0]*root[1]*root[2]);
    b = std::real(root[0]*root[1]+root[0]*root[2]+root[1]*root[2]);
    c = -std::real(root[0]+root[1]+root[2]);
}
void coeffs_from_roots(double& a, double&b, double&c,
                       const std::vector<double>& real_root)
{
    assert(real_root.size()==3);
    std::array<std::complex<double>, 3> root;
    for(size_t i=0; i<root.size(); i++) {
        root[i] = real_root[i];
    }
    coeffs_from_roots(a,b,c,root);
}

// random coefficients for a cubic polynomial: f(x) = a + b*x + c*x^2 + x^3
// returns vector of all distinct real roots if known
// input type:
//  - 1: three real roots
//  - 2: one double, one single real root
//  - 3: one triple real root
//  - 4: one real and two conjugate complex root
//  - 5: random coefficients (without knowing the exact roots)
std::vector<double> rand_coeffs(int type, double& a, double&b, double&c)
{
    std::array<std::complex<double>, 3> complex_root;
    if(type==1) {
        // three real roots
        for(size_t i=0; i<complex_root.size(); i++) {
            complex_root[i] = rand_laplace();
        }
    } else if(type==2) {
        // one double, one single real root
        complex_root[0] = rand_laplace();
        complex_root[1] = rand_laplace();
        complex_root[2] = complex_root[1];
    } else if(type==3) {
        // one triple real root
        complex_root[0] = rand_laplace();
        complex_root[1] = complex_root[0];
        complex_root[2] = complex_root[0];
    } else if(type==4) {
        // one real and two conjugate complex root
        complex_root[0] = std::complex<double>(rand_laplace(),rand_laplace());
        complex_root[1] = std::conj(complex_root[0]);
        complex_root[2] = rand_laplace();
    } else if(type==5) {
        // random coefficients (without knowing the exact roots)
        a = rand_laplace();
        b = rand_laplace();
        c = rand_laplace();
    }

    std::vector<double> root;
    if(type!=5) {
        // calculate polynomial coeffs
        coeffs_from_roots(a,b,c,complex_root);

        // return real roots only
        for(size_t i=0; i<complex_root.size(); i++) {
            if(std::imag(complex_root[i])==0.0) {
                root.push_back(std::real(complex_root[i]));
            }
        }
    }
    // sort and erase duplicates
    std::sort(root.begin(),root.end());
    root.erase(std::unique(root.begin(), root.end()), root.end());

    return root;
}
// global variable to describe the type of roots
std::vector<std::string> type_name = {
    "", "three roots",
    "two roots (one double)", "one root (triple)",
    "one root, two complex roots",  "random coefficients"
};



// helper function to do min/max and sum for errors in: x and f(x)
void collect_stats(const std::vector<double>& exact,
                   const std::vector<double>& numeric,
                   double a, double b, double c,
                   double& sum_xerr, double& max_xerr,
                   double& sum_yerr, double& max_yerr,
                   int& count, int& count_incorrect)
{
    if(exact.size()>0) {
        if(numeric.size()!=exact.size()) {
            count_incorrect++;
        } else {
            for(size_t j=0; j<numeric.size(); j++) {
                double xerr = fabs(numeric[j]-exact[j]);
                double yerr = fabs(f(a,b,c,numeric[j]));
                sum_xerr += xerr;
                sum_yerr += yerr;
                count++;
                if(xerr>max_xerr)
                    max_xerr=xerr;
                if(yerr>max_yerr)
                    max_yerr=yerr;
            }
        }
    } else {
        // no information about exact roots
        for(size_t j=0; j<numeric.size(); j++) {
            double yerr = fabs(f(a,b,c,numeric[j]));
            sum_yerr += yerr;
            count++;
            if(yerr>max_yerr)
                max_yerr=yerr;
        }
    }
}

// test correctness of results using num_loops different examples per case
// - randomly generate coefficients based on known or unknown roots
//   - 4 different cases where roots are known and coeffs calculated
//   - 1 case where coefficients are randomly generated (roots unknown)
// - compare against known roots and check the residual error f(x)
void run_comparison(int num_loops)
{
    srand(0);
    double a,b,c;                       // polynomial coefficients
    std::vector<double> x;              // numerical roots
    std::vector<double> exact;          // exact roots

    printf("           x-error        f(x)-error    wrong no of roots\n");
    printf("          avg    max      avg    max\n");
    for(int type=1; type<=5; type++) {
        // count cases when number of roots wrong (e.g. one missing or extra)
        std::array<int,3> count= {0}, count_incorrect= {0};
        // average and max error in the solution
        std::array<double,3> avg_xerr= {0}, max_xerr= {0};
        // average and max error in the function value
        std::array<double,3> avg_yerr= {0}, max_yerr= {0};

        for(int i=0; i<num_loops; i++) {
            exact = rand_coeffs(type,a,b,c);
            x = gsl_poly_wrapper(a,b,c,1.0);
            collect_stats(exact,x,a,b,c, avg_xerr[0],max_xerr[0],
                          avg_yerr[0],max_yerr[0], count[0],count_incorrect[0]);

            x = tk::internal::solve_cubic(a,b,c,1.0,0);
            collect_stats(exact,x,a,b,c, avg_xerr[1],max_xerr[1],
                          avg_yerr[1],max_yerr[1], count[1],count_incorrect[1]);

            x = tk::internal::solve_cubic(a,b,c,1.0,1);
            collect_stats(exact,x,a,b,c, avg_xerr[2],max_xerr[2],
                          avg_yerr[2],max_yerr[2], count[2],count_incorrect[2]);
        }
        std::array<std::string,3> name = {"gsl  ", "tk(0)", "tk(1)"};
        printf("type %i: %s\n", type, type_name[type].c_str());
        for(size_t i=0; i<avg_xerr.size(); i++) {
            avg_xerr[i]/=(double)count[i];
            avg_yerr[i]/=(double)count[i];
            printf("%s:  %6.0e %6.0e,  %6.0e %6.0e, %10i\n",
                   name[i].c_str(),
                   avg_xerr[i], max_xerr[i], avg_yerr[i], max_yerr[i],
                   count_incorrect[i]);
        }
    }
}



// benchmark run for the gsl cubic root finding routine
double bench_gsl(int num_loops, int type, int newton_iter=0)
{
    (void)newton_iter;      // unused variable
    double a,b,c;
    double x0=0.0, x1=0.0, x2=0.0;
    srand(0);

    double sum=0.0;
    for(int i=0; i<num_loops; i++) {
        if(i%100==0) {
            rand_coeffs(type,a,b,c);        // new coefficients
        }
        gsl_poly_solve_cubic(c,b,a,&x0,&x1,&x2);
        a+=1e-100;   // just so the loop is not optimised away
        sum+=x0;
    }
    return sum;
}

// benchmark run for tk::internal::solve_cubic()
double bench_tk(int num_loops, int type, int newton_iter=0)
{
    double a,b,c;
    srand(0);

    double sum=0.0;
    for(int i=0; i<num_loops; i++) {
        if(i%100==0) {
            rand_coeffs(type,a,b,c);        // new coefficients
        }
        std::vector<double> x=tk::internal::solve_cubic(a,b,c,1.0,newton_iter);
        a+=1e-100;   // just so the loop is not optimised away
        sum+=x[0];
    }
    return sum;
}


// cpu clock frequency
double mhz;

#define RUN_BENCH(DESCR,ROUTINE,OPS,TYPE,ITER)              \
{                                                           \
   printf("%-35s loops=%.0e,", (DESCR).c_str(), (double) OPS);           \
      double t=stoptime();                                  \
      double y=ROUTINE(OPS,TYPE,ITER);                      \
      t=stoptime()-t;                                       \
      printf("\t%.3fs ",t);                                 \
      double cycles=t/(OPS)*mhz*1e6;                        \
      if(cycles<10000)  printf("(%4.0f cycl)", cycles);     \
      else              printf("(%.1e cycl)", cycles);      \
      fflush(stdout);                                       \
      if(y==-123.456789123) printf("only to make use of y");  \
   printf("\n");                                            \
}

void run_bench(double loops)
{
    for(int type=1; type<=5; type++)
        RUN_BENCH(std::string("gsl:   ")+type_name[type], bench_gsl, loops, type, 0);
    printf("\n");
    for(int type=1; type<=5; type++)
        RUN_BENCH(std::string("tk:    ")+type_name[type], bench_tk, loops, type, 0);
    printf("\nwith one newton iteration:\n");
    for(int type=1; type<=5; type++)
        RUN_BENCH(std::string("tk(1): ")+type_name[type], bench_tk, loops, type, 1);
}


void usage(char** argv)
{
    printf("root finding for f(x)= a + b*x + c*x^2 + d*x^3: \n");
    printf("usage: \n");
    printf(" a) %s -coeffs <a> <b> <c> <d>\n", argv[0]);
    printf(" b) %s -roots <x0> <x1> <x2>\n", argv[0]);
    printf(" c) %s -bench <cpu mhz> <loops in thousands>\n", argv[0]);
    printf(" d) %s -test <loops in thousands>\n", argv[0]);
    exit(EXIT_FAILURE);
}

int main(int argc, char** argv)
{

    if(argc<3 || argc>6) {
        usage(argv);
    }

    std::string param = argv[1];
    if(param=="-coeffs" || param=="-roots") {
        // calculate roots for a single cubic plynomial
        double a,b,c,d;
        std::array<std::vector<double>,5> sol;  // solutions x to f(x)=0
        std::array<std::string,5> name;         // solution method name
        if(param=="-coeffs") {
            if(argc!=6)
                usage(argv);
            a=atof(argv[2]);
            b=atof(argv[3]);
            c=atof(argv[4]);
            d=atof(argv[5]);
        } else { // param=="-roots"
            if(argc!=5)
                usage(argv);
            sol[0].resize(3);
            sol[0][0]=atof(argv[2]);            // 1st real root
            sol[0][1]=atof(argv[3]);            // 2nd real root
            sol[0][2]=atof(argv[4]);            // 3rd real root
            coeffs_from_roots(a,b,c,sol[0]);
            d=1.0;
            // sort and remove duplicates
            std::sort(sol[0].begin(),sol[0].end());
            sol[0].erase(std::unique(sol[0].begin(),sol[0].end()),sol[0].end());
        }
        printf("f(x) = %.16f + %.16f*x + %.16f*x^2 + %.16f*x^3\n",a,b,c,d);
        name[0]="exact";
        name[1]="gsl";
        sol[1] = gsl_poly_wrapper(a, b, c, d);
        name[2]="tk(0)";
        sol[2] = tk::internal::solve_cubic(a, b, c, d, 0);  // no newton iter
        name[3]="tk(1)";
        sol[3] = tk::internal::solve_cubic(a, b, c, d, 1);  // 1 newton iter
        name[4]="tk(2)";
        sol[4] = tk::internal::solve_cubic(a, b, c, d, 2);  // 2 newton iter
        printf("\nsolutions x to f(x)=0\n");
        for(size_t i=0; i<sol.size(); i++) {
            printf("%-5s: ", name[i].c_str());
            for(size_t j=0; j<sol[i].size(); j++) {
                printf("%20.16f ", sol[i][j]);
            }
            printf("\n");
        }
        printf("\nf(x)\n");
        for(size_t i=0; i<sol.size(); i++) {
            printf("%-5s: ", name[i].c_str());
            for(size_t j=0; j<sol[i].size(); j++) {
                printf("%20.2e ", f(a,b,c,d,sol[i][j]));
            }
            printf("\n");
        }

    } else if(param=="-bench") {
        // calculate roots for many cubic plynomials and measure time
        if(argc!=4)
            usage(argv);
        mhz=atof(argv[2]);
        int loops=atoi(argv[3])*1e3;
        run_bench(loops);

    } else if(param=="-test") {
        // calculate roots for many cubic plynomials and compare accuracy
        int loops=atoi(argv[2])*1e3;
        run_comparison(loops);

    } else {
        printf("unknown parameter: %s\n", param.c_str());
        usage(argv);
    }

    printf("\n");
    return EXIT_SUCCESS;
}


