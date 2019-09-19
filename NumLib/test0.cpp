#include <iostream>

#include "Numerics/Arg/ArgN.h"
#include "Numerics/Ivp/Utils.h"
#include "Numerics/Equ/Utils.h"


int main()
{
	using Vec = Num::Arg::VecN<double, 2>;
	using Mat = Num::Arg::MatNxM<double, 2, 2>;

	auto func = [&] (double t, double u)
	{
		return t * u * u;
	};

	auto jacob = Num::Ivp::make_diff_derivative<double, double>(func, 1e-5);

	auto f = [&] (double t)
	{
		return std::exp(t);
	};

	auto deriv = Num::Equ::make_diff_derivative(f, 1e-6);

	std::cout << jacob(1.0, 5.0) << std::endl;
	std::cout << deriv(2.0) << std::endl;

	std::cin.get();
}