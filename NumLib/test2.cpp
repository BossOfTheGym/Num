#include <iostream>
#include <array>
#include <cmath>
#include <string>
#include <algorithm>

#include "Utility/Utility.h"
#include "Elliptic/Elliptic.h"
#include "Parabolic/Parabolic.h"


using Rect = Utility::Rect<double>;

using Template2x3 = Utility::Template<double, 2, 3>;
using Template3x3 = Utility::Template<double, 3, 3>;

using EllipticSolverSimple     = Elliptic::SolverRect3x3_1_Simple<double>;
using EllipticSolverZeidel     = Elliptic::SolverRect3x3_1_Zeidel<double>;
using EllipticSolverRelaxation = Elliptic::SolverRect3x3_1_Relaxation<double>;

using ParabolicSolver = Parabolic::SolverRect2x3_1<double>;


template<class Solution, class TestFunction>
double testSolution(int nt, int nx, const Solution& numerical, TestFunction&& solution)
{
	double diff = 0.0;

	for (int i = 0; i <= nt; i++)
	{
		for (int j = 0; j <= nx; j++)
		{
			if (std::isnan(numerical[i][j]))
			{
				std::cout << "NaN" << std::endl;
			}

			diff = std::max(diff, std::abs(numerical[i][j] - solution(i, j)));
		}
	}

	return diff;
}


template<class Solution, class Stream>
void writeSolutionToStream(
	const Solution& solution
	, int nt, int nx
	, double x0, double t0
	, double ht, double hx
	, Stream& os
)
{
	os << nt << " " << nx << std::endl;
	os << t0 << " " << x0 << std::endl;
	os << ht << " " << hx << std::endl;

	for (auto& elem : solution)
	{
		os << elem << " ";
	}
}


template<class Solution, class Stream>
void writeSolutionToStreamCompressed(
	const Solution& solution
	, int nt, int nx
	, double x0, double t0
	, double ht, double hx
	, int ratio_t, int ratio_x
	, Stream& os
)
{
	os << nt / ratio_t << " " << nx / ratio_x << std::endl;

	os << t0 << " " << x0 << std::endl;

	os << ht * ratio_t << " " << hx * ratio_x << std::endl;

	for (int i = 0; i <= nt; i += ratio_t)
	{
		for (int j = 0; j <= nx; j += ratio_x)
		{
			os << solution[i][j] << " ";
		}
	}
}



void testElliptic()
{
	const double x0 = 0.0;
	const double x1 = 1.0;
	const double t0 = 0.0;
	const double t1 = 1.0;

	const double dx = x1 - x0;
	const double dt = t1 - t0;

	const int nx = 100;
	const int nt = 100;

	const double hx = dx / nx;
	const double ht = dt / nt;

	const double hx2 = hx * hx;
	const double ht2 = ht * ht;

	
	//basic functions
	auto function = [&] (double t, double x)
	{
		return 2.0 * std::exp(-(t + x));
	};

	
	auto leftBoundFunction = [&] (double t)
	{
		return -std::exp(-t);
	};

	auto rightBoundFunction = [&] (double t)
	{
		return std::exp(-(t + 1.0));
	};

	auto lowerBoundFunction = [&] (double x)
	{
		return std::exp(-x);
	};

	auto upperBoundFunction = [&] (double x)
	{
		return std::exp(-(1.0 + x));
	};


	auto leftLowerFunction = [&] ()
	{
		return 1.0;
	};

	auto rightLowerFunction = [&] ()
	{
		return std::exp(-1.0);
	};

	auto leftUpperFunction = [&] ()
	{
		return std::exp(-1.0);
	};

	auto rightUpperFunction = [&] ()
	{
		return std::exp(-2.0);
	};


	//equations
	auto equation2 = [&] (int i, int j)
	{
		Template3x3 equ;

		double a = 1.0 / ht2;
		double b = 1.0 / hx2;
		double c = 0.0;

		equ[2][0] = c;     equ[2][1] = a - 2.0 * c;        equ[2][2] = c;
		equ[1][0] = b + c; equ[1][1] = -2.0 * (a + b + c); equ[1][2] = b + c;
		equ[0][0] = c;     equ[0][1] = a - 2.0 * c;        equ[0][2] = c;

		return equ;
	};

	auto equation4 = [&] (int i, int j)
	{
		Template3x3 equ;

		double a = 1.0 / ht2;
		double b = 1.0 / hx2;
		double c = (a + b) / 12.0;

		equ[2][0] = c;         equ[2][1] = a - 2 * c;            equ[2][2] = c;
		equ[1][0] = b - 2 * c; equ[1][1] = -2 * (a + b) + 4 * c; equ[1][2] = b - 2 * c;
		equ[0][0] = c;         equ[0][1] = a - 2 * c;            equ[0][2] = c;

		return equ;
	};


	auto leftBound2 = [&] (int i)
	{
		Template3x3 equ;

		double a = 1.0 / hx;
		double b = hx / 2 / ht2;

		equ[2][0] = 0.0; equ[2][1] = b;          equ[2][2] = 0.0;
		equ[1][0] = 0.0; equ[1][1] = -a - 2 * b; equ[1][2] = a;
		equ[0][0] = 0.0; equ[0][1] = b;          equ[0][2] = 0.0;

		return equ;
	};

	auto leftBound4 = [&] (int i)
	{
		Template3x3 equ;

		double a = 1 / hx;
		double b = hx / (2 * ht2);
		double c = (ht2 + hx2) / ht2 / hx / 12;

		equ[2][0] = 0.0; equ[2][1] = b - c; equ[2][2] = c;
		equ[1][0] = 0.0; equ[1][1] = -a - 2 * b + 2 * c; equ[1][2] = a - 2 * c;
		equ[0][0] = 0.0; equ[0][1] = b - c; equ[0][2] = c;

		return equ;
	};


	auto rightBound = [&] (int i)
	{
		Template3x3 equ;

		equ[2][0] = 0.0; equ[2][1] = 0.0; equ[2][2] = 0.0;
		equ[1][0] = 0.0; equ[1][1] = 1.0; equ[1][2] = 0.0;
		equ[0][0] = 0.0; equ[0][1] = 0.0; equ[0][2] = 0.0;

		return equ;
	};

	auto lowerBound = [&] (int j)
	{
		Template3x3 equ;

		equ[2][0] = 0.0; equ[2][1] = 0.0; equ[2][2] = 0.0;
		equ[1][0] = 0.0; equ[1][1] = 1.0; equ[1][2] = 0.0;
		equ[0][0] = 0.0; equ[0][1] = 0.0; equ[0][2] = 0.0;

		return equ;
	};

	auto upperBound = [&] (int j)
	{
		Template3x3 equ;

		equ[2][0] = 0.0; equ[2][1] = 0.0; equ[2][2] = 0.0;
		equ[1][0] = 0.0; equ[1][1] = 1.0; equ[1][2] = 0.0;
		equ[0][0] = 0.0; equ[0][1] = 0.0; equ[0][2] = 0.0;

		return equ;
	};


	auto leftLowerCorner = [&] ()
	{
		Template3x3 equ;

		equ[2][0] = 0.0; equ[2][1] = 0.0; equ[2][2] = 0.0;
		equ[1][0] = 0.0; equ[1][1] = 1.0; equ[1][2] = 0.0;
		equ[0][0] = 0.0; equ[0][1] = 0.0; equ[0][2] = 0.0;

		return equ;
	};

	auto rightLowerCorner = [&] ()
	{
		Template3x3 equ;

		equ[2][0] = 0.0; equ[2][1] = 0.0; equ[2][2] = 0.0;
		equ[1][0] = 0.0; equ[1][1] = 1.0; equ[1][2] = 0.0;
		equ[0][0] = 0.0; equ[0][1] = 0.0; equ[0][2] = 0.0;

		return equ;
	};

	auto leftUpperCorner = [&] ()
	{
		Template3x3 equ;

		equ[2][0] = 0.0; equ[2][1] = 0.0; equ[2][2] = 0.0;
		equ[1][0] = 0.0; equ[1][1] = 1.0; equ[1][2] = 0.0;
		equ[0][0] = 0.0; equ[0][1] = 0.0; equ[0][2] = 0.0;

		return equ;
	};

	auto rightUpperCorner = [&] ()
	{
		Template3x3 equ;

		equ[2][0] = 0.0; equ[2][1] = 0.0; equ[2][2] = 0.0;
		equ[1][0] = 0.0; equ[1][1] = 1.0; equ[1][2] = 0.0;
		equ[0][0] = 0.0; equ[0][1] = 0.0; equ[0][2] = 0.0;

		return equ;
	};


	//approximations of Higher Order
	auto functionHO2 = [&] (int i, int j)
	{
		double t = t0 + ht * i;
		double x = x0 + hx * j;

		return function(t, x);
	};

	auto functionHO4 = [&] (int i, int j)
	{
		double f12 = function(t0 + ht * i, x0 + hx * (j + 1));
		double f10 = function(t0 + ht * i, x0 + hx * (j - 1));

		double f11 = function(t0 + ht * i, x0 + hx * j);

		double f21 = function(t0 + ht * (i + 1), x0 + hx * j);
		double f01 = function(t0 + ht * (i - 1), x0 + hx * j);

		return 2.0 * f11 / 3.0 + (f21 + f01 + f12 + f10) / 12.0;
	};

	
	auto leftBoundFunctionHO2 = [&] (int i)
	{
		double t = t0 + ht * i;

		return leftBoundFunction(t) + hx / 2 * function(t, x0);
	};

	auto leftBoundFunctionHO4 = [&] (int i)
	{
		double lb2 = leftBoundFunction(t0 + ht * (i + 1));
		double lb1 = leftBoundFunction(t0 + ht * (i    ));
		double lb0 = leftBoundFunction(t0 + ht * (i - 1));


		double f12 = function(t0 + ht * i, x0 + hx * (+1));
		double f10 = function(t0 + ht * i, x0 + hx * (-1));

		double f11 = function(t0 + ht * i, x0);

		double f21 = function(t0 + ht * (i + 1), x0);
		double f01 = function(t0 + ht * (i - 1), x0);

		return lb1 
			+ (ht2 - hx2) / (6.0 * ht2) * (lb2 - 2.0 * lb1 + lb0)
			+ hx / 2.0 * (2.0 / 3.0 * f11 + f21 / 12.0 + f01 / 12.0 - f10 / 12.0 + f12 / 4.0);
	};


	auto rightBoundFunctionHO = [&] (int i)
	{
		double t = t0 + ht * i;

		return rightBoundFunction(t);
	};

	auto lowerBoundFunctionHO = [&] (int j)
	{
		double x = x0 + hx * j;

		return lowerBoundFunction(x);
	};

	auto upperBoundFunctionHO = [&] (int j)
	{
		double x = x0 + hx * j;

		return upperBoundFunction(x);
	};


	auto leftLowerFunctionHO = [&] ()
	{
		return leftLowerFunction();
	};

	auto rightLowerFunctionHO = [&] ()
	{
		return rightLowerFunction();
	};

	auto leftUpperFunctionHO = [&] ()
	{
		return leftUpperFunction();
	};

	auto rightUpperFunctionHO = [&] ()
	{
		return rightUpperFunction();
	};


	//exact solution
	auto solution = [&] (double t, double x)
	{
		return std::exp(-(t + x));
	};

	auto solutionNet = [&] (int i, int j)
	{
		double t = t0 + ht * i;
		double x = x0 + hx * j;

		return solution(t, x);
	};


	//solver lambda + cache
	auto cache = Rect(nt + 1, nx + 1, 0.0);

	auto solve2 = [&] (auto&& solver) -> decltype(auto)
	{
		const int limit = nt * nx * 2;
		const double eps = ht2 * hx2 / 100;

		return Elliptic::solveCached(
			cache, solver

			, nt, nx, eps, limit

			, equation2, functionHO2
			, leftBound2, leftBoundFunctionHO2
			, rightBound, rightBoundFunctionHO
			, lowerBound, lowerBoundFunctionHO
			, upperBound, upperBoundFunctionHO

			, rightLowerCorner, rightLowerFunctionHO
			, rightUpperCorner, rightUpperFunctionHO
			, leftLowerCorner, leftLowerFunctionHO
			, leftUpperCorner, leftUpperFunctionHO
		);
	};

	auto solve4 = [&] (auto&& solver) -> decltype(auto)
	{
		const int limit = nt * nx * 2;
		const double eps = ht2 * hx2 / 100;

		return Elliptic::solveCached(
			cache, solver

			, nt, nx, eps, limit

			, equation4, functionHO4
			, leftBound4, leftBoundFunctionHO4
			, rightBound, rightBoundFunctionHO
			, lowerBound, lowerBoundFunctionHO
			, upperBound, upperBoundFunctionHO

			, rightLowerCorner, rightLowerFunctionHO
			, rightUpperCorner, rightUpperFunctionHO
			, leftLowerCorner, leftLowerFunctionHO
			, leftUpperCorner, leftUpperFunctionHO
		);
	};


	//solvers
	//4-th order
	//Jacoby iterations
	auto solver1 = EllipticSolverSimple();
	auto numerical = solve4(solver1);

	auto diff1 = testSolution(nt, nx, numerical, solutionNet);

	//Gauss-Zeidel method
	auto solver2 = EllipticSolverZeidel();
	numerical = solve4(solver2);

	auto diff2 = testSolution(nt, nx, numerical, solutionNet);

	//relaxation
	auto solver3 = EllipticSolverRelaxation(1.5);
	numerical = solve4(solver3);

	auto diff3 = testSolution(nt, nx, numerical, solutionNet);

	//result
	std::cout << "4-th order errors" << std::endl;
	std::cout << diff1 << " " << diff2 << " " << diff3 << std::endl;
	std::cin.get();

	//2-nd order
	//Jacoby iterations
	numerical = solve2(solver1);

	diff1 = testSolution(nt, nx, numerical, solutionNet);

	//Gauss-Zeidel method	
	numerical = solve2(solver2);

	diff2 = testSolution(nt, nx, numerical, solutionNet);

	//relaxation
	numerical = solve2(solver3);

	diff3 = testSolution(nt, nx, numerical, solutionNet);

	//result
	std::cout << "2-th order errors" << std::endl;
	std::cout << diff1 << " " << diff2 << " " << diff3 << std::endl;
	std::cin.get();
}


void testParabolic()
{	
	const double PI  = std::acos(-1);
	const double SQ2 = std::sqrt(2);

	const double t0 = 0.0;
	const double t1 = 1.0;
	const double x0 = 0.0;
	const double x1 = 1.0;
	
	const double dt = t1 - t0;
	const double dx = x1 - x0;
	
	const int nt = 10000;
	const int nx = 10000;
	
	const double ht = dt / nt;
	const double hx = dx / nx;
	
	const double hx2 = hx * hx;

	//basic functions
	//TEST
	auto function = [&] (double t, double x)
	{
		return SQ2 * std::sin(PI / 4 + t + x);
	};

	//TEST
	auto leftBoundFunction = [&] (double t)
	{
		return std::sin(t);
	};


	auto rightBoundFunction = [&] (double t)
	{
		return std::sin(t + 1.0);
	};

	auto lowerBoundFunction = [&] (double x)
	{
		return std::sin(x);
	};


	//equations
	//TEST
	auto equation2 = [&] (int i, int j)
	{
		Template2x3 equ;

		double alpha = 0.0;
		double sigma = 0.5 + alpha * hx2 / ht;

		double a = sigma / hx2;
		double b = (1.0 - sigma) / hx2;
		double c = 1.0 / ht;

		equ[1][0] = -a; equ[1][1] = +c + 2 * a; equ[1][2] = -a;
		equ[0][0] = -b; equ[0][1] = -c + 2 * b; equ[0][2] = -b;

		return equ;
	};

	//TEST
	auto equation4 = [&] (int i, int j)
	{
		Template2x3 equ;

		double sigma = 0.5 - hx2 / (12.0 * ht);

		double a = sigma / hx2;
		double b = (1.0 - sigma) / hx2;
		double c = 1.0 / ht;

		equ[1][0] = -a; equ[1][1] = +c + 2 * a; equ[1][2] = -a;
		equ[0][0] = -b; equ[0][1] = -c + 2 * b; equ[0][2] = -b;

		return equ;
	};


	auto leftBound2 = [&] (int i)
	{
		Template2x3 equ;

		equ[1][0] = 0.0; equ[1][1] = 1.0; equ[1][2] = 0.0;
		equ[0][0] = 0.0; equ[0][1] = 0.0; equ[0][2] = 0.0;

		return equ;
	};

	auto leftBound4 = [&] (int i)
	{
		Template2x3 equ;

		equ[1][0] = 0.0; equ[1][1] = 1.0; equ[1][2] = 0.0;
		equ[0][0] = 0.0; equ[0][1] = 0.0; equ[0][2] = 0.0;

		return equ;
	};


	auto rightBound = [&] (int i)
	{
		Template2x3 equ;

		equ[1][0] = 0.0; equ[1][1] = 1.0; equ[1][2] = 0.0;
		equ[0][0] = 0.0; equ[0][1] = 0.0; equ[0][2] = 0.0;

		return equ;
	};


	//approximations of Higher Order
	auto functionHO2 = [&] (int i, int j)
	{
		double f21 = function(t0 + ht * (i + 1), x0 + hx * j);
		double f11 = function(t0 + ht * (i    ), x0 + hx * j);
		double f01 = function(t0 + ht * (i - 1), x0 + hx * j);

		return f11 + (f21 - f01) / 4;
	};
	auto functionHO4 = [&] (int i, int j)
	{
		double f21 = function(t0 + ht * (i + 1), x0 + hx * j);

		double f12 = function(t0 + ht * i, x0 + hx * (j + 1));
		double f11 = function(t0 + ht * i, x0 + hx * (j    ));
		double f10 = function(t0 + ht * i, x0 + hx * (j - 1));

		return f10 / 12 + f11 / 3 + f12 / 12 + f21 / 2;
	};

	
	auto leftBoundFunctionHO2 = [&] (int i)
	{
		return leftBoundFunction(t0 + ht * i);
	};
	
	auto leftBoundFunctionHO4 = [&] (int i)
	{
		return leftBoundFunction(t0 + ht * i);
	};


	auto rightBoundFunctionHO = [&] (int i)
	{
		return rightBoundFunction(t0 + ht * i);
	};

	auto lowerBoundFunctionHO = [&] (int j)
	{
		return lowerBoundFunction(x0 + hx * j);
	};


	//exact solution
	auto solution = [&] (double t, double x)
	{
		return std::sin(t + x);
	};

	auto solutionNet = [&] (int i, int j)
	{
		return solution(t0 + ht * i, x0 + hx * j);
	};

	//solver lambda
	auto solve2 = [&] (auto&& solver)
	{
		return solver.solve(
			nt, nx
			, equation2 , functionHO2
			, leftBound2, leftBoundFunctionHO2
			, rightBound, rightBoundFunctionHO
			            , lowerBoundFunctionHO
		);
	};

	auto solve4 = [&] (auto&& solver)
	{
		return solver.solve(
			nt, nx
			, equation4 , functionHO4
			, leftBound4, leftBoundFunctionHO4
			, rightBound, rightBoundFunctionHO
			            , lowerBoundFunctionHO
		);
	};

	//test
	auto solver = ParabolicSolver();

	auto numerical = solve4(solver);
	auto diff1 = testSolution(nt, nx, numerical, solutionNet);
}


void testBalance()
{

}

double test()
{
	std::vector<double> a(10000000);

	double b = 0.0;
	for (int k = 0; k < 10; k++)
	{
		for (int i = 0; i < 10000000; i++)
		{
			a[i] = 1.0;
			b += a[i];
			a[i] += b;
		}

	}

	auto p = &a;
	auto pb = &b;

	std::cout << p << std::endl;
	std::cout << pb << std::endl;

	return *pb;
}

int main()
{
	testElliptic();
	
	//testParabolic();

	return 0;
}