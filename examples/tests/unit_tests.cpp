#include "spline.h"

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <vector>
#include <string>


#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SplineUnitTests
#include <boost/test/unit_test.hpp>


// setup different types of splines for testing
// which are at least "cont_deriv" times continuously differentiable
std::vector<tk::spline> setup_splines(const std::vector<double>& X,
                                      const std::vector<double> Y,
                                      int cont_deriv,
                                      bool make_all_monotonic = false)
{
    assert(X.size()==Y.size());
    std::vector<tk::spline> splines;

    if(cont_deriv<=2) {
        // C^2 splines
        // default cubic spline
        splines.push_back(tk::spline());
        splines.back().set_points(X,Y);
        // default cubic spline: non-zero curvature at boundary
        splines.push_back(tk::spline());
        splines.back().set_boundary(tk::spline::second_deriv, 0.2,
                                    tk::spline::second_deriv, -0.3);
        splines.back().set_points(X,Y);
    }
    if(cont_deriv<=1) {
        // C^1 splines
        // default cubic spline: left derivative given
        splines.push_back(tk::spline());
        splines.back().set_boundary(tk::spline::first_deriv, 0.5,
                                    tk::spline::second_deriv, 0.0);
        splines.back().set_points(X,Y);
        // default cubic spline: right derivative given
        splines.push_back(tk::spline());
        splines.back().set_boundary(tk::spline::second_deriv, 0.0,
                                    tk::spline::first_deriv, 1.2);
        splines.back().set_points(X,Y);
        // default cubic spline: left and right derivative given
        splines.push_back(tk::spline());
        splines.back().set_boundary(tk::spline::first_deriv, 0.0,
                                    tk::spline::first_deriv, -1.2);
        splines.back().set_points(X,Y);

        // cubic hermite spline
        splines.push_back(tk::spline());
        splines.back().set_points(X,Y,tk::spline::cspline_hermite);
        // cubic hermite spline: left derivative given
        splines.push_back(tk::spline());
        splines.back().set_boundary(tk::spline::first_deriv, 0.5,
                                    tk::spline::second_deriv, 0.0);
        splines.back().set_points(X,Y,tk::spline::cspline_hermite);
        // cubic hermite spline: right derivative given
        splines.push_back(tk::spline());
        splines.back().set_boundary(tk::spline::second_deriv, 0.0,
                                    tk::spline::first_deriv, 1.2);
        splines.back().set_points(X,Y,tk::spline::cspline_hermite);
        // cubic hermite spline: left and right derivative given
        splines.push_back(tk::spline());
        splines.back().set_boundary(tk::spline::first_deriv, 0.0,
                                    tk::spline::first_deriv, -1.2);
        splines.back().set_points(X,Y,tk::spline::cspline_hermite);
    }
    if(cont_deriv<=0) {
        // C^0 splines
        // linear interpolation
        splines.push_back(tk::spline());
        splines.back().set_points(X,Y,tk::spline::linear);
    }

    // add the same splines but with monotinicity requested
    if(cont_deriv<=1 && make_all_monotonic==false) {
        size_t n=splines.size();
        for(size_t i=0; i<n; i++) {
            splines.push_back(splines[i]);
            splines.back().make_monotonic();
        }
    }
    // make all monotonic
    if(make_all_monotonic==true) {
        size_t n=splines.size();
        for(size_t i=0; i<n; i++) {
            splines[i].make_monotonic();
        }
    }

    return splines;
}


std::vector<double> evaluation_grid(const std::vector<double>& X,
                                    double margin, size_t more_points,
                                    bool include_vector=true)
{
    std::vector<double> eval_grid;
    if(include_vector==true) {
        eval_grid=X;
    }
    for(size_t i=0; i<more_points; i++) {
        double p = X[0]-margin + (X.back()-X[0]+2.0*margin)*i/more_points;
        eval_grid.push_back(p);
    }
    return eval_grid;
}


// checks spline is C^0, i.e. continuous and goes through all grid points
BOOST_AUTO_TEST_CASE( Continuity )
{
    const double dx = 1e-15;        // step size to check continuity
    const double max_deriv = 10.0;  // max 1st deriv of below spline
    // using well posed grid data for spline interpolation
    std::vector<double> X = {-10.2,  0.1,  0.8, 2.4,  3.1,  3.2, 4.7, 19.1};
    std::vector<double> Y = {  2.7, -1.2, -0.5, 1.5, -1.0, -0.7,  1.2, -1.3};
    // points to check continuity of spline
    std::vector<double> eval_grid = evaluation_grid(X,1.0,1000);

    // setup all possible types of splines which are at least C^0
    std::vector<tk::spline> splines=setup_splines(X,Y,0);

    for(size_t k=0; k<splines.size(); k++) {
        const tk::spline& s=splines[k];
        BOOST_TEST_CONTEXT("spline " << k << ": " << s.info())
        {
        for(size_t i=0; i<eval_grid.size(); i++) {
            double x = eval_grid[i];
            if(i<Y.size()) {
                double y = Y[i];
                BOOST_CHECK_CLOSE( s(x), y, 0.0 );          // hit grid points
            }
            BOOST_CHECK_SMALL(s(x)-s(x+dx), max_deriv*dx);  // left continuous
            BOOST_CHECK_SMALL(s(x)-s(x-dx), max_deriv*dx);  // right continous
        }
        } // BOOST_TEST_CONTEXT
    }
}

// checks spline is C^1, i.e. differentiable
BOOST_AUTO_TEST_CASE( Differentiability )
{
    const double dx = 1e-15;        // step size to check continuity
    const double h  = 2e-8;         // step size for finite differences
    const double max_deriv = 100.0; // max 2nd deriv of below spline
    std::vector<double> X = {-10.2,  0.1,  0.8, 2.4,  3.1,  3.2, 4.7, 19.1};
    std::vector<double> Y = {  2.7, -1.2, -0.5, 1.5, -1.0, -0.7,  1.2, -1.3};
    // points to check continuity/differentiability of spline
    std::vector<double> eval_grid = evaluation_grid(X,1.0,1000);

    // setup all possible types of splines which are at least C^1
    std::vector<tk::spline> splines=setup_splines(X,Y,1);

    for(size_t k=0; k<splines.size(); k++) {
        const tk::spline& s=splines[k];
        BOOST_TEST_CONTEXT("spline " << k << ": " << s.info())
        {
        for(size_t i=0; i<eval_grid.size(); i++) {
            double x = eval_grid[i];
            double right_diff =   (s(x+h)-s(x)) / h;
            double left_diff =    (s(x)-s(x-h)) / h;
            double central_diff = (s(x+h)-s(x-h) ) / (2.0*h);
    
            BOOST_CHECK_SMALL( left_diff-right_diff, max_deriv*h );
            BOOST_CHECK_SMALL( central_diff-right_diff, max_deriv*h);
            BOOST_CHECK_SMALL( s.deriv(1,x)-central_diff, max_deriv*h );
    
            BOOST_CHECK_SMALL( s.deriv(1,x-dx)-s.deriv(1,x), max_deriv*dx );
            BOOST_CHECK_SMALL( s.deriv(1,x+dx)-s.deriv(1,x), max_deriv*dx );
        }
        } // BOOST_TEST_CONTEXT
    }
}


// checks spline is C^2
BOOST_AUTO_TEST_CASE( Smoothness )
{
    const double dx = 1e-15;        // step size to check continuity
    const double h  = 3e-6;         // step size for finite differences
    const double max_deriv = 1000.0;    // max 3nd deriv of below spline
    std::vector<double> X = {-10.2,  0.1,  0.8, 2.4,  3.1,  3.2, 4.7, 19.1};
    std::vector<double> Y = {  2.7, -1.2, -0.5, 1.5, -1.0, -0.7,  1.2, -1.3};
    // points to check continuity/differentiability of spline
    std::vector<double> eval_grid = evaluation_grid(X,1.0,1000);

    // setup all possible types of splines which are at least C^2
    std::vector<tk::spline> splines=setup_splines(X,Y,2);

    for(size_t k=0; k<splines.size(); k++) {
        const tk::spline& s=splines[k];
        BOOST_TEST_CONTEXT("spline " << k << ": " << s.info())
        {
        for(size_t i=0; i<eval_grid.size(); i++) {
            double x = eval_grid[i];
            double right_diff =   (s(x)-2.0*s(x+h)+s(x+2.0*h)) / (h*h);
            double left_diff =    (s(x-2.0*h)-2.0*s(x-h)+s(x)) / (h*h);
            double central_diff = (s(x-h)-2.0*s(x)+s(x+h)) / (h*h);
    
            BOOST_CHECK_SMALL( left_diff-right_diff, max_deriv*h );
            BOOST_CHECK_SMALL( central_diff-right_diff, max_deriv*h );
            BOOST_CHECK_SMALL( s.deriv(2,x)-central_diff, max_deriv*h );
    
            BOOST_CHECK_SMALL( s.deriv(2,x-dx)-s.deriv(2,x), max_deriv*dx );
            BOOST_CHECK_SMALL( s.deriv(2,x+dx)-s.deriv(2,x), max_deriv*dx );
        }
        } // BOOST_TEST_CONTEXT
    }
}

// checks the calculation of the third derivative
BOOST_AUTO_TEST_CASE( ThirdDeriv )
{
    // 2nd derivative is piecewise linear, so we can use a larger
    // step size for finite differences for 3rd derivative
    const double h    = 1e-5;     // step size for finite differences
    const double tol  = 1e-10;    // error tolerance in 3rd derivative
    std::vector<double> X = {-10.2,  0.1,  0.8, 2.4,  3.1,  3.2, 4.7, 19.1};
    std::vector<double> Y = {  2.7, -1.2, -0.5, 1.5, -1.0, -0.7,  1.2, -1.3};

    // evaluation grid must not include any of the grid points X
    std::vector<double> eval_grid = { -12.3, -4.3, 0.5, 3.5, 7.3, 22.1, 95.0};

    // setup all possible types of splines which are at least C^0
    std::vector<tk::spline> splines=setup_splines(X,Y,0);

    for(size_t k=0; k<splines.size(); k++) {
        const tk::spline& s=splines[k];
        BOOST_TEST_CONTEXT("spline " << k << ": " << s.info())
        {
        for(size_t i=0; i<eval_grid.size(); i++) {
            double x = eval_grid[i];
            double right_diff =   (s.deriv(2,x+h)-s.deriv(2,x)) / h;
            double left_diff =    (s.deriv(2,x)-s.deriv(2,x-h)) / h;
            double central_diff = (s.deriv(2,x+h)-s.deriv(2,x-h) ) / (2.0*h);
    
            BOOST_CHECK_SMALL( left_diff-right_diff, tol );
            BOOST_CHECK_SMALL( central_diff-right_diff, tol );
            BOOST_CHECK_SMALL( s.deriv(3,x)-central_diff, tol );
    
            BOOST_CHECK_SMALL( s.deriv(3,x-h)-s.deriv(3,x), 0.0 );  // piecew. const
            BOOST_CHECK_SMALL( s.deriv(3,x+h)-s.deriv(3,x), 0.0 );  // piecew. const
        }
        } // BOOST_TEST_CONTEXT
    }
}

BOOST_AUTO_TEST_CASE( Monotonicity )
{
    std::vector<double> X = {-10.2,  0.1, 0.8,  2.4,  3.1,  3.2,  4.7, 19.1};
    std::vector<double> Y = {   1.0, 1.0, 0.0, -0.1, -0.2, -1.0, -1.1, -1.2};
    // points to check monotinicity of the spline
    std::vector<double> eval_grid = evaluation_grid(X,0.0,2000,false);

    // setup all possible types of splines which are at least C^1 and monotonic
    std::vector<tk::spline> splines=setup_splines(X,Y,1,true);

    for(size_t k=0; k<splines.size(); k++) {
        const tk::spline& s=splines[k];
        BOOST_TEST_CONTEXT("spline " << k << ": " << s.info())
        {
        int monotonic=0;
        for(size_t i=1; i<eval_grid.size(); i++) {
            double x1=eval_grid[i-1];
            double x2=eval_grid[i];
            if(monotonic==0) {
                // determin whether decreasing or increasing
                if(s(x1)<s(x2)) {
                    monotonic=1;
                } else if(s(x1)>s(x2)) {
                    monotonic=-1;
                }
            } else if(monotonic==-1) {
                // decreasing
                BOOST_CHECK( s(x1)>=s(x2) );
            } else if(monotonic==1){
                // increasing
                BOOST_CHECK( s(x1)<=s(x2) );
            } else {
                assert(false);
            }
        }
        } // BOOST_TEST_CONTEXT
    }
}



BOOST_AUTO_TEST_CASE( BoundaryTypes )
{
    const double dx = 1e-5;        // step size to check continuity
    std::vector<double> X = {-10.2,  0.1,  0.8, 2.4,  3.1,  3.2, 4.7, 19.1};
    std::vector<double> Y = {  2.7, -1.2, -0.5, 1.5, -1.0, -0.7,  1.2, -1.3};
    {
        // natural boundary conditions, i.e. 2nd derivative is zero
        tk::spline s;
        s.set_points(X,Y);
        BOOST_CHECK_CLOSE( s.deriv(2,X[0]-dx), 0.0, 0.0 );      // 2nd deriv = 0
        BOOST_CHECK_CLOSE( s.deriv(2,X.back()+dx), 0.0, 0.0 );  // 2nd deriv = 0
    }

    {
        // left and right with given 2nd and 1st derivative, respectively
        tk::spline s;
        double deriv1=1.3;
        double deriv2=-0.7;
        s.set_boundary(tk::spline::second_deriv, deriv2,
                       tk::spline::first_deriv, deriv1);
        s.set_points(X,Y);
        BOOST_CHECK_CLOSE( s.deriv(2,X[0]), deriv2, 0.0 );
        BOOST_CHECK_CLOSE( s.deriv(2,X[0]-1.0), deriv2, 0.0 );
        BOOST_CHECK_CLOSE( s.deriv(1,X.back()), deriv1, 1e-12);
        BOOST_CHECK_CLOSE( s.deriv(1,X.back()+1.0), deriv1, 1e-12);
        BOOST_CHECK_CLOSE( s.deriv(2,X.back()+1.0), 0.0, 0.0);
    }
    {
        // left and right with given 1st and 2nd derivative, respectively
        tk::spline s;
        double deriv1=-4.7;
        double deriv2=-1.2;
        s.set_boundary(tk::spline::first_deriv, deriv1,
                       tk::spline::second_deriv, deriv2);
        s.set_points(X,Y);
        BOOST_CHECK_CLOSE( s.deriv(2,X.back()), deriv2, 0.0 );
        BOOST_CHECK_CLOSE( s.deriv(2,X.back()+1.0), deriv2, 0.0 );
        BOOST_CHECK_CLOSE( s.deriv(1,X[0]), deriv1, 1e-12);
        BOOST_CHECK_CLOSE( s.deriv(1,X[0]-1.0), deriv1, 1e-12);
        BOOST_CHECK_CLOSE( s.deriv(2,X[0]-1.0), 0.0, 0.0);
    }
}


// function to be approximated by a spline
double myfunc(double x)
{
    return sin(x);
}

// this test checks how well a smooth function is approximated by a spline
BOOST_AUTO_TEST_CASE( FunctionApproximation )
{
    double lo=0.0;
    double hi=2.0*M_PI;
    // regression results
    std::vector< std::vector<double> > max_errors = {
        {
            1.0, 1.087710e-01, 2.001701e-02, 6.866929e-04,
            3.189426e-05, 1.764343e-06, 1.043467e-07, 6.352541e-09,
            3.919548e-10, 2.432621e-11
        },
        {
            1.0, 1.895384e-01, 9.289323e-02, 7.816602e-03,
            1.082736e-03, 1.319855e-04, 1.607416e-05, 1.977196e-06,
            2.450258e-07, 3.049079e-08
        }
    };
    std::vector< std::vector<double> > avg_errors = {
        {
            6.365561e-01, 4.767218e-02, 1.161861e-02, 2.479389e-04,
            1.100872e-05, 6.015310e-07, 3.546300e-08, 2.157237e-09,
            1.330842e-10, 8.264881e-12
        },
        {
            6.365561e-01, 9.277315e-02, 3.244986e-02, 3.617337e-03,
            3.095662e-04, 3.116561e-05, 3.514112e-06, 4.189284e-07,
            5.122710e-08, 6.337061e-09
        }
    };
    double tol=1e-6;        // 6 correct digits

    // evaluation grid
    std::vector<double> eval_grid(10000);
    for(size_t i=0; i<eval_grid.size(); i++)
        eval_grid[i] = lo + (double)i/(eval_grid.size()-1) * (hi-lo);

    // spline sizes
    std::vector<size_t> N = {3, 4, 5, 10, 20, 40, 80, 160, 320, 640};

    // spline types
    std::vector<tk::spline::spline_type> types =
    {tk::spline::cspline, tk::spline::cspline_hermite};
    std::vector<std::string> type_names = {"cspline", "hermite"};
    for(size_t l=0; l<2; l++) {
        tk::spline::spline_type type=types[l];
        const std::vector<double>& avg_error = avg_errors[l];
        const std::vector<double>& max_error = max_errors[l];
        const std::string type_name = type_names[l];

        for(size_t k=0; k<N.size(); k++) {
            size_t n=N[k];
            std::vector<double> X(n), Y(n);
            for(size_t i=0; i<n; i++) {
                X[i] = lo + (double)i/(n-1) * (hi-lo);
                Y[i] = myfunc(X[i]);
            }
            tk::spline s;
            s.set_points(X,Y,type);
            for(size_t i=0; i<n; i++) {
                X[i] = lo + (double)i/(n-1) * (hi-lo);
                Y[i] = myfunc(X[i]);
            }
            double avg = 0.0;
            double max = 0.0;
            for(size_t i=0; i<eval_grid.size(); i++) {
                double x = eval_grid[i];
                double err = std::fabs(myfunc(x)-s(x));
                avg += err;
                if(err>max)
                    max=err;
            }
            avg /= eval_grid.size();
            printf("approx %s: n = %3lu, avg error = %e, max error = %e\n",
                   type_name.c_str(), n, avg, max );
            BOOST_CHECK_CLOSE(avg, avg_error[k], tol*100.0);
            BOOST_CHECK_CLOSE(max, max_error[k], tol*100.0);
        }
    }
}

double rand_unif()
{
    return (double)std::rand()/RAND_MAX;
}

// this test is intended to detect any numerical changes
// and will require re-basing whenever a numerical change is expected
BOOST_AUTO_TEST_CASE( RegressionTest )
{
    std::vector<size_t> N = {3, 4, 5, 9, 13, 22, 86, 137, 657, 12357};

    double sum=0.0;
    double sum_mod=0.0;
    for(size_t k=0; k<N.size(); k++) {
        size_t n=N[k];
        std::vector<double> X(n), Y(n);
        X[0]=rand_unif();
        Y[0]=rand_unif();
        for(size_t i=1; i<n; i++) {
            X[i] = X[i-1] + rand_unif();            // monotonic increasing
            Y[i] = Y[i-1] + rand_unif()-0.5;
        }

        // evaluation grid
        std::vector<double> eval_grid(1000);
        for(size_t i=0; i<eval_grid.size(); i++)
            eval_grid[i] = X[0]-5.0 + rand_unif()*(X.back()-X[0]+10.0);
        eval_grid.insert(eval_grid.end(), X.begin(), X.end());  // add X points

        // setup all possible types of splines which are at least C^0
        std::vector<tk::spline> splines=setup_splines(X,Y,0);

        for(size_t l=0; l<splines.size(); l++) {
            const tk::spline& s = splines[l];
            for(size_t i=0; i<eval_grid.size(); i++) {
                double x = eval_grid[i];
                double val = s(x)+0.1*s.deriv(1,x)
                             + 0.01*s.deriv(2,x)+0.001*s.deriv(3,x);
                sum += val;
                sum_mod += std::fmod(val, 1.0);
                sum_mod = std::fmod(sum_mod, 1.0);
            }
            //printf("spline %lu size %5lu: sum=%.16e, mod=%.16f\n",l,n,sum,sum_mod);
        }
    }
    //printf("sum=%.16e, mod=%.16f\n", sum, sum_mod);
    BOOST_CHECK_CLOSE(sum, -298638197797.79303, 0.0);
    BOOST_CHECK_CLOSE(sum_mod, 0.20285718270622555, 0.0);
}
