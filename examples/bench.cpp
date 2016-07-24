#include <cstdio>
#include <cstdlib>
#include <time.h>
#include <sys/time.h>
#include "spline.h"           // tk
#include "interpolation.h"    // alglib



#define RUN_BENCH(DESCR,ROUTINE,OPS)                        \
{                                                           \
   printf("%s loops=%.0e,", DESCR, (double) OPS);           \
   for(int i=bench::TK;  i<=bench::ALGLIB; i++){            \
      bench::method m = static_cast<bench::method>(i);      \
      double t=stoptime();                                  \
      double y=ROUTINE(m,OPS);                              \
      t=stoptime()-t;                                       \
      printf("\t%.3fs ",t);                                 \
      double cycles=t/(OPS)*mhz*1e6;                        \
      if(cycles<10000)  printf("(%4.0f cycl)", cycles);     \
      else              printf("(%.1e cycl)", cycles);      \
      fflush(stdout);                                       \
      if(y==-123.456789123) printf("only to make use of y");  \
   }                                                        \
   printf("\n");                                            \
}

double stoptime(void)
{
    struct timeval t;
    gettimeofday(&t,NULL);
    return (double) t.tv_sec + t.tv_usec/1000000.0;
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





// gnu libc has a reasonable rand() implementation
double rand01()
{
    return rand()/(RAND_MAX+1.0);
}


class bench
{
private:
    alglib::spline1dinterpolant m_spline_al;
    tk::spline m_spline_tk;
    std::vector<double> m_X1, m_Y1, m_X2, m_Y2;
    alglib::real_1d_array m_AX1, m_AY1, m_AX2, m_AY2;

public:
    enum method {TK, ALGLIB};
    // init vectors: X1 and X2 random grid points, in approximately [0,1]
    void init_vectors(size_t n1, size_t n2)
    {
        const double dx1=1.0/(0.5*n1);
        const double dx2=(dx1*n1)/n2;
        m_X1.resize(n1);
        m_Y1.resize(n1);
        m_X1[0]=0.0;
        m_Y1[0]=rand01();
        for(size_t i=1; i<n1; i++) {
            m_X1[i]=m_X1[i-1]+rand01()*dx1;
            m_Y1[i]=rand01();
        }
        m_X2.resize(n2);
        m_X2[0]=0.0;
        for(size_t i=1; i<n2; i++) {
            m_X2[i]=m_X2[i-1]+rand01()*dx2;
        }

        // copy data to alglib arrays
        m_AX1.setcontent(n1, &(m_X1[0]));
        m_AY1.setcontent(n1, &(m_Y1[0]));
        m_AX2.setcontent(n2, &(m_X2[0]));
    }

    // create spline from (X,Y) vectors
    void create_spline(method m)
    {
        if(m==TK) {
            m_spline_tk.set_points(m_X1,m_Y1);
        } else {
            alglib::spline1dbuildcubic(m_AX1, m_AY1, m_X1.size(),2,0.0,2,0.0,
                                       m_spline_al);
        }
    }

    // spline(x)
    double get_value(method m, double x)
    {
        double res;
        if(m==TK) {
            res=m_spline_tk(x);
        } else {
            res=alglib::spline1dcalc(m_spline_al,x);
        }
        return res;
    }

    // convert values from one grid to another (X1,Y1) --> (X2,?)
    void convert(method m)
    {
        if(m==TK) {
            tk::spline s;
            size_t n2=m_X2.size();
            s.set_points(m_X1,m_Y1);
            m_Y2.resize(n2);
            for(size_t i=0; i<n2; i++) {
                m_Y2[i]=s(m_X2[i]);
            }
        } else {
            spline1dconvcubic(m_AX1, m_AY1, m_X1.size(),2,0.0,2,0.0,
                              m_AX2, m_X2.size(), m_AY2);
        }
    }

    // benchmark spline creation
    double bench_create_spline(method m, int num_ops)
    {
        for(int i=0; i<num_ops; i++) {
            create_spline(m);
        }
        return get_value(m,0.0);
    }
    // benchmark spline evaluation at a specifc point x: spline(x)
    double bench_get_value(method m, int num_ops)
    {
        create_spline(m);
        double x=0.5;
        for(int i=0; i<num_ops; i++) {
            double y=get_value(m,x);
            // generate a fairly random sequence in [0,1]
            x+=y;
            x-=floor(x);
        }
        return x;
    }
    // benchmark grid conversion (X1,X2) --> (X2, ?)
    double bench_convert(method m, int num_ops)
    {
        for(int i=0; i<num_ops; i++) {
            convert(m);
        }
        size_t n=m_X2.size()-1;
        return m==TK ? m_Y2[n] : m_AY2[n];
    }

    // compare accuracy between two methods of spline(x)
    void compare(double& l2, double& linf, method m1, method m2, int num_ops)
    {
        double a=m_X1[0];
        double b=m_X1.back();
        l2=0.0;
        linf=0.0;
        for(int i=0; i<num_ops; i++) {
            double x=a+(b-a)*rand01();
            double y1=get_value(m1,x);
            double y2=get_value(m2,x);
            double err=fabs(y1-y2);
            l2+=err*err;
            linf=std::max(linf,err);
        }
        l2=sqrt(l2)/num_ops;
    }

};



// cpu clock frequency
double mhz;


// ---------------------------------------------------------------------
// main
// ---------------------------------------------------------------------
int main(int argc, char** argv)
{
    if(argc<4) {
        printf("usage: %s <cpu mhz> <spline size> <loops in thousand>\n",
               argv[0]);
        exit(EXIT_FAILURE);
    }
    mhz=atof(argv[1]);
    int dim=atoi(argv[2]);
    int loops=atoi(argv[3])*1e3;
    int loops2 = (int) round_dec((double)loops/dim,0);

    bench b;
    b.init_vectors(dim,dim);

    printf("\t\t\t\ttk\t\t\talglib\n");
    RUN_BENCH("random access:  ", b.bench_get_value,loops);
    RUN_BENCH("spline creation:", b.bench_create_spline,loops2);
    RUN_BENCH("grid transform: ", b.bench_convert, loops2);

    double l2, linf;
    b.compare(l2,linf,bench::TK, bench::ALGLIB, loops);
    printf("\naccuracy: max difference = %.2e, l2-norm difference = %.2e\n",
           linf, l2);

    return EXIT_SUCCESS;
}


