#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <functional>


typedef std::function<double(double)> Function;


double intSimpson_1_3(const Function& func, double a, double b, int n)
{
    const double h = (b - a) / n;

    double integral = 0.0;

    for (int i = 0; i < n; i++)
    {
        integral += 4 * func(a + h * (i + 0.5));
    }
    for (int i = 1; i < n; i++)
    {
        integral += 2 * func(a + h * i);
    }

    return (integral + func(a) + func(b)) * h / 6;
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

double rungeRule1(
      double int1, double h1
    , double int2, double h2
    , double order)
{
    return (int2 - int1) / (1.0 - pow(h2 / h1, order));
}

double rungeRule2(
      double int1, double h1
    , double int2, double h2
    , double order)
{
    return (int2 - int1) / (pow(h1 / h2, order) - 1.0);
}


double intGaussRule_2(
	const Function_1& func
	, double a
	, double b
	, int n
)
{
	//Gauss args
	const double X0 = -1.0 / sqrt(3);
	const double X1 = +1.0 / sqrt(3);

	const double h = (b - a) / n;

	double integral = 0.0;

	const double A = h / 2;
	for (int i = 0; i < n; i++)
	{
		const double B = a + h * (i + 0.5);

		integral += func(A * X0 + B);
		integral += func(A * X1 + B);
	}

	return integral * A;
}


double intGaussRule_2(
	const Function_2& func
	, double a, double b
	, double c, double d
	, int n
)
{
	//Gauss args
	const double X0 = -1.0 / sqrt(3);
	const double X1 = +1.0 / sqrt(3);

	const double h0 = (b - a) / n;
	const double h1 = (d - c) / n;

	const double A0 = h0 / 2;
	const double A1 = h1 / 2;

	double integral = 0.0;

	for (int i = 0; i < n; i++)
	{
		const double B0 = a + h0 * (i + 0.5);

		for (int j = 0; j < n; j++)
		{
			const double B1 = c + h1 * (j + 0.5);

			integral += func(A0 * X0 + B0, A1 * X0 + B1);
			integral += func(A0 * X0 + B0, A1 * X1 + B1);
			integral += func(A0 * X1 + B0, A1 * X0 + B1);
			integral += func(A0 * X1 + B0, A1 * X1 + B1);
		}
	}

	return integral * A0 * A1;
}

int main()
{
    auto func = [] (double x) -> double
    {
        return log10(sqrt(x * x + 2) + 1) 
            / (cos(x * x - 3 * x + 5) + 2);
    };

    double a = -1.0;
    double b = +5.0;

    int n1 = 2500;
    int n2 = 5000;

    double intSimpson1 = intSimpson_1_3(func, -1.0, 5.0, n1 / 2);
    double intSimpson2 = intSimpson_1_3(func, -1.0, 5.0, n2 / 2);

    double intTrapezoidal1 = intTrapezoidalRule(func, -1.0, 5.0, n1);
    double intTrapezoidal2 = intTrapezoidalRule(func, -1.0, 5.0, n2);

    double h1 = (b - a) / 2500;
    double h2 = (b - a) / 5000;

    double rungeSimpson     = rungeRule2(intSimpson1, 2 * h1, intSimpson2, 2 * h2, 4);
    double rungeTrapezoidal = rungeRule2(intTrapezoidal1, h1, intTrapezoidal2, h2, 1);

    std::cout.precision(15);

    std::cout << "Simpson: "   << intSimpson2  << std::endl;
    std::cout << "Runge rule:" << rungeSimpson << std::endl << std::endl;

    std::cout << "Trapezoidal: " << intTrapezoidal2  << std::endl;
    std::cout << "Runge rule:"   << rungeTrapezoidal << std::endl << std::endl;

    std::cin.get();
    std::cin.get();

    return 0;
}