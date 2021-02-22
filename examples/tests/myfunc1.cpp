// test that we can include spline.h in more than one compilation unit
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "spline.h"

class myclass
{
protected:
    tk::spline m_spline1;
    tk::spline m_spline2;
public:
    myclass(const std::vector<double>& X, const std::vector<double>& Y):
        m_spline1(X,Y)
    {
        std::vector<double> X2 = {0.0, 0.5, 1.0, 2.0};
        std::vector<double> Y2 = {0.6, 0.5, 0.3, 1.1};
        m_spline2.set_points(X2,Y2);
    }
    double getval(double x) const
    {
        return m_spline1(x)+m_spline2(x);
    }
};

double myfunc1(double x)
{

    static std::vector<double> X = {0.1, 0.4, 1.2, 1.8, 2.0};
    static std::vector<double> Y = {0.1, 0.7, 0.6, 1.1, 0.9};
    static myclass s(X,Y);
    return s.getval(x);
}
