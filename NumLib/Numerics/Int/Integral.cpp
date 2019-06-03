#include <cmath>
#include <functional>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <iomanip>

typedef std::function<double(double)> Function;

const double PI_2 = asin(1);
const double PI   = 2 * PI_2;


//lazy version
double intSimpson_1_3(const Function& func, double a, double b, int n)
{
    const double h = (b - a) / n;

    double integral = 0.0;

    for (int i = 0; i < n; i++)
    {
        const double val0 = func(a + h * i);
        const double val1 = func(a + h * (i + 0.5));
        const double val2 = func(a + h * (i + 1.0));

        integral += val0 + 4 * val1 + val2;
    }

    return integral * h / 6;
}
//lazy version
double intSimpson_3_8(const Function& func, double a, double b, int n)
{
    const double h = (b - a) / n;

    double integral = 0.0;

    for (int i = 0; i < n; i++)
    {
        const double val0 = func(a + h * i);
        const double val1 = func(a + h * (i + 1.0 / 3));
        const double val2 = func(a + h * (i + 2.0 / 3));
        const double val3 = func(a + h * (i + 1.0));

        integral += val0 + 3 * val1 + 3 * val2 + val3;
    }

    return integral * h / 8;
}


double intLeftRectangles(const Function& func, double a, double b, int n)
{
    const double h = (b - a) / n;

    double integral = 0.0;

    for (int i = 0; i < n; i++)
    {
        const double val0 = func(a + h * i);

        integral += val0;
    }

    return integral * h;
}

double intRightRectangles(const Function& func, double a, double b, int n)
{
    const double h = (b - a) / n;

    double integral = 0.0;

    for (int i = 0; i < n; i++)
    {
        const double val0 = func(a + h * (i + 1));

        integral += val0;
    }

    return integral * h;
}

double intMidpointRectangles(const Function& func, double a, double b, int n)
{
    const double h = (b - a) / n;

    double integral = 0.0;

    for (int i = 0; i < n; i++)
    {
        const double val0 = func(a + h * (i + 1.0 / 2));

        integral += val0;
    }

    return integral * h;
}


double intTrapezoidalRule(const Function& func, double a, double b, int n)
{
    const double h = (b - a) / n;

    double integral = 0.0;

    for (int i = 1; i < n; i++)
    {
        const double val0 = func(a + h * i);

        integral += val0;
    }

    return (integral + (func(a) + func(b)) / 2) * h;
}


int main()
{
    auto func = [] (double x) -> double
    {
        return x / sqrt(x * x + 4 * x + 1);
    };

    const double a = 1;
    const double b = 2;

    const int n0 = 5;
    const int n1 = 7;

    std::cout.precision(12);

    std::cout << std::left << std::showpos;

    std::cout << std::setw(20) << "Left rule" << ": " << intLeftRectangles(func, a, b, n0) << std::endl
        << std::setw(20) << "Right rule"     << ": " << intRightRectangles(func, a, b, n0)    << std::endl
        << std::setw(20) << "Midpoint rule"  << ": " << intMidpointRectangles(func, a, b, n0) << std::endl

        << std::setw(20) << "Trapezoidal rule" << ": " << intTrapezoidalRule(func, a, b, n0) << std::endl

        << std::setw(20) << "1/3 rule" << ": " << intSimpson_1_3(func, a, b, n0) << std::endl
        << std::setw(20) << "3/8 rule" << ": " << intSimpson_3_8(func, a, b, n1) << std::endl;

    std::cout << std::endl << "Computing finished" << std::endl;

    std::cin.get();

    return 0;
}