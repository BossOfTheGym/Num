#include <iostream>
#include <sstream>

#include "Numerics/Arg/ArgN.h"
#include "Numerics/Equ/Utils.h"
#include "Numerics/Ivp/Utils.h"
#include "Numerics/Ivp/RungeKutta.h"
#include "Numerics/Ivp/ImplicitSolverAdaptor.h"


void testIvpUtils()
{
	std::ostringstream output;

	using Vec = Num::Arg::VecN<double, 2>;
	using Mat = Num::Arg::MatNxM<double, 2, 2>;

	auto func = [&] (double t, double u)
	{
		return t * u * u;
	};

	auto jacob = Num::Ivp::make_diff_derivative<double>(func, 1e-5);

	auto f = [&] (double t)
	{
		return std::exp(t);
	};

	auto deriv = Num::Equ::make_diff_derivative(f, 1e-6);

	output << "---Ivp utils test---" << std::endl;
	output << jacob(1.0, 5.0)[0][0] << std::endl;
	output << deriv(2.0)[0][0] << std::endl;
	output << "---Ivp utils test finished" << std::endl << std::endl;

	std::cout << output.str();
}

void testRke()
{
	using Vec = Num::Arg::VecN<double, 1>;
	using Mat = Num::Arg::MatNxM<double, 1, 1>;

	auto func = [&] (double t, double u)
	{
		return u;
	};

	auto solver = Num::Ivp::make_rke_solver<1, double, double>(Num::Ivp::ExplicitMethods<double>::classic4());

	double h = 0.0001;
	double t0 = 0.0;
	double up = 1.0;
	double un;

	for (int i = 0; i < 10000; i++)
	{
		un = solver.solve(func, t0 + i * h, up, h);
		up = un;
	}
	std::cout << t0 + 10000 * h << std::endl;
	std::cout << "Num: " << up << " Exact: " << std::exp(1.0) << std::endl;
}

void testRki()
{
	using Vec = Num::Arg::VecN<double, 1>;
	using Mat = Num::Arg::MatNxM<double, 1, 1>;

	auto func = [&] (double t, Vec u)
	{
		return Vec(-u);
	};

	auto solver = Num::Ivp::make_rki_solver<1, double, Vec>(
		Num::Ivp::ImplicitMethods<double>::midpoint2(), 1e-12, 20
	);
	auto adapter = Num::Ivp::make_implicit_adaptor<double, Vec>(solver);

	double h = 0.0001;
	double t0 = 0.0;
	Vec up = 1.0;
	Vec un;

	for (int i = 0; i < 10000; i++)
	{
		un = adapter.solve(func, t0 + i * h, up, h);
		up = un;
	}
	std::cout << t0 + 10000 * h << std::endl;
	std::cout << "Num: " << up[0] << " Exact: " << std::exp(-1.0) << std::endl;
}


int main()
{
	testRke();
	testRki();
	std::cin.get();
}