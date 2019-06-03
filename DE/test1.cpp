#include <vector>
#include <iostream>
#include <functional>
#include <vector>


#include "Numerics/Arg/ArgN.h"
#include "Numerics/LinAlg/Gauss.h"
#include "Numerics/Utils/Utils.h"
#include "Numerics/Equ/Neuton.h"
#include "Numerics/Equ/FixedPoint.h"
#include "Numerics/Ivp/RungeKutta.h"

/*#include <iostream>
#include <vector>
#include <utility>
#include <algorithm>

#include "Numerics/Arg/ArgN.h"
#include "Numerics/Ivp/RungeKutta.h"



using Num::Ivp::RungeKuttaExplicit;
using Num::Ivp::RungeKuttaImplicit;
using Num::Ivp::Methods;

using Mat = Num::Arg::MatNxM<double, 3, 3>;
using Y   = Num::Arg::VecN<double, 3>;


struct Node
{
	double t;
	Y y;
};

template<int N>
struct Helper
{
	static constexpr int value = N;
};


template<class Function, class ExplicitSolver>
auto solveWithExplicit(ExplicitSolver&& solver, Function&& ode, double t0, Y y0, double h, int count)
{
	auto solution = std::vector<Node>(count + 1);

	solution[0] = Node{t0, y0};
	for (int i = 1; i <= count; i++)
	{
		auto lastY = solution[i - 1].y;
		auto lastT = (t0 + (i - 1) * h);

		auto[t, y] = solver.solve(ode, lastT, lastY, h);

		solution[i] = Node{(t0 + i * h), y};
	}

	return solution;
}

template<class Function, class Jacobian, class ImplicitSolver>
auto solveWithImplicit(ImplicitSolver&& solver, Function&& ode, Jacobian&& jac, double t0, Y y0, double h, int count)
{
	auto solution = std::vector<Node>(count + 1);

	solution[0] = Node{t0, y0};
	for (int i = 1; i <= count; i++)
	{
		auto lastY = solution[i - 1].y;
		auto lastT = (t0 + (i - 1) * h);

		auto[t, y] = solver.solve(ode, jac, lastT, lastY, h);

		solution[i] = Node{(t0 + i * h), y};
	}

	return solution;
}


template<class Tableau, class Function, class Order>
auto solveWithRKE(Tableau&& tableau, Function&& ode, double t0, Y y0, double h, int count, Order&& dummy)
{
	auto solver = RungeKuttaExplicit<Order::value, double, Y, decltype(std::move(tableau))>(std::move(tableau));

	return solveWithExplicit(solver, ode, t0, y0, h, count);
}

template<class Function, class Jacobian, class Tableau, class Order>
auto solveWithRKI(Tableau&& tableau, Function&& ode, Jacobian&& jac, double t0, Y y0, double h, int count, Order&& dummy)
{
	auto solver = RungeKuttaImplicit<Order::value, double, Y, decltype(std::move(tableau))>(std::move(tableau), 1e-16, 50);

	return solveWithImplicit(solver, ode, jac, t0, y0, h, count);
}


template<class Function, class Jacobian, class Order>
void test(Function&& ode, Jacobian&& jac, double t0, Y y0, double h, int count, Order&& dummy)
{
	//explicit
	auto classic4 = Methods<double>::classic4();
	auto ralston4 = Methods<double>::ralston4();
	auto rule38_4 = Methods<double>::three_eighths_rule4();

	//implicit
	auto backwardsEuler = Methods<double>::backwordEuler1();
	auto midpoint2      = Methods<double>::midpoint2();
	auto gaussLegendre4 = Methods<double>::gaussLegendre4();
	auto gaussLegendre6 = Methods<double>::gaussLegendre6();


	std::vector<Node> solution;

	solution = solveWithRKE(std::move(classic4), std::move(ode), t0, y0, h, count, std::move(dummy));
	solution = solveWithRKE(std::move(ralston4), std::move(ode), t0, y0, h, count, std::move(dummy));;
	solution = solveWithRKE(std::move(rule38_4), std::move(ode), t0, y0, h, count, std::move(dummy));

	solution = solveWithRKI(std::move(backwardsEuler), std::move(ode), std::move(jac), t0, y0, h, count, std::move(dummy));
	solution = solveWithRKI(std::move(midpoint2)     , std::move(ode), std::move(jac), t0, y0, h, count, std::move(dummy));
	solution = solveWithRKI(std::move(gaussLegendre4), std::move(ode), std::move(jac), t0, y0, h, count, std::move(dummy));
	solution = solveWithRKI(std::move(gaussLegendre6), std::move(ode), std::move(jac), t0, y0, h, count, std::move(dummy));
}

void firstOde()
{
	auto t0 = double(0.0);
	auto y0 = Y(0.05, 0.05, 0.9);
	auto h = 1e6;
	auto count = 1000000;

	auto ode = [] (double t, Y y)
	{
		double y1 = y[0];
		double y2 = y[1];
		double y3 = y[2];

		return Y(
			-0.04 * y1 + 1e4 * y2 * y3
			, +0.04 * y1 - 1e4 * y2 * y3 - 3e7 * y2 * y2
			, +3e7 * y2 * y2
		);
	};

	auto jac = [] (double t, Y y)
	{
		double y1 = y[0];
		double y2 = y[1];
		double y3 = y[2];

		return Mat{
			-0.04, +1e4 * y3           , +1e4 * y2
			, +0.04, -1e4 * y3 - 6e7 * y2, -1e4 * y2
			, +0.0 ,            +6e7 * y2, +0.0
		};
	};

	auto invariant = [] (Y y)
	{
		double y1 = y[0];
		double y2 = y[1];
		double y3 = y[2];

		return y1 + y2 + y3;
	};

	test(std::move(ode), std::move(jac), t0, y0, h, count, Helper<3>());

	std::cin.get();
}


void secondOde()
{

}

int main()
{
	firstOde();

	return 0;
}*/

void testLinAlg()
{

    using Num::Arg::MatNxM;
    using Num::Arg::VecN;
    using Num::LinAlg::GaussEliminationSingle;

    MatNxM<double, 3, 3> mat{
          1.0, 5.0, 7.0
        , 0.0, 2.0, 3.0
        , 1.0, 1.0, 0.0
    };

    VecN<double, 3> solution{1.0, 2.0, 3.0};

    VecN<double, 3> term = mat * solution;

    GaussEliminationSingle<double, 3> solver;

    solver.solve(mat, term);

    std::cout << "*** Testing LinAlg module ***" << std::endl;
    std::cout 
        << "Solution expected: " 
        << solution[0] << ", " 
        << solution[1] << ", " 
        << solution[2] << std::endl;
    std::cout << "Solution got: (" 
        << term[0] << ", " 
        << term[1] << ", " 
        << term[2] << " )" << std::endl;
    std::cout << "*** LinAlg test finished ***" << std::endl << std::endl;
}

void testNeuton()
{
    using Num::EPS;
    using Num::Equ::NeutonScalar;

    auto neutonSolver = NeutonScalar<double>(10, 1e-12);
    auto root = neutonSolver.solve(
        [] (double x) -> double
        {
            return x * x - 9.0;
        }
        ,[] (double x) -> double
        {
            return 2 * x;
        }
        , 0.5
    );

    using Num::Arg::VecN;
    using Num::Arg::MatNxM;
    using Num::Utils::Norm;
    using Num::Utils::NormType;
    using Num::Equ::NeutonSystem;
    using UsedNorm = Norm<VecN<double, 2>, NormType::SPHERE>;
	using Num::Equ::DiffJacobian;

	using Vec = VecN<double, 2>;
	using Mat = MatNxM<double, 2, 2>;

    NeutonSystem<double, 2> systemSolver(10, EPS<double>);

	auto func = [] (const Vec& arg) -> Vec
	{
		const double& x = arg[0];
		const double& y = arg[1];

		return VecN<double, 2>{x * x + y * y - 1.0, x - y};
	};
    auto rootSys = systemSolver.solve(
        std::move(func)
        , DiffJacobian<Vec, Mat, decltype(func)>(std::move(func), 1e-10)
        , Vec{1.0, 0.0}
    );


    std::cout << "Neuton solver test" << std::endl;
    std::cout << "Equation root expected: " << 3 << std::endl;
    std::cout << "Got       : " << root.first << std::endl;
    std::cout << "Iterations: " << root.second << std::endl << std::endl;

    std::cout << "Equation system root expected: (" << 1.0 / sqrt(2.0) << ", " << 1.0 / sqrt(2.0) << ")" << std::endl;
    std::cout << "Got       : (" << rootSys.first[0] << ", " << rootSys.first[1] << ")" << std::endl;
    std::cout << "Iterations: " << rootSys.second << std::endl;
    std::cout << "Neuton solver test finished" << std::endl << std::endl;
}

void testFixed()
{
    using Num::EPS;
    using Num::Equ::FixedPointScalar;
    
    auto root = FixedPointScalar<double>(20).solve(
        [] (double x) -> double
        {
            //x^2 - 2 = 0
            return x - (x * x - 2.0) / 4.0;
        }
        , 1.5
    );


    using Num::Arg::VecN;
    using Num::Arg::MatNxM;
    using Num::Equ::FixedPointSystem;

    auto rootSys = FixedPointSystem<double, 2>(20).solve(
        [] (const VecN<double, 2>& argument) -> VecN<double, 2>
        {
            const double& x = argument[0];
            const double& y = argument[1];

            return VecN<double, 2>({y, sqrt(x + 2.0)});
        }
        , VecN<double, 2>({0.0, 1.0})
    );


    std::cout << "FixedPoint solver test" << std::endl;
    std::cout << "Equation root expected: " << sqrt(2) << std::endl;
    std::cout << "Got       : " << root.first << std::endl;
    std::cout << "Iterations: " << root.second << std::endl << std::endl;

    std::cout << "Equation system root expected: (" << 2 << ", " << 2 << ")" << std::endl;
    std::cout << "Got       : (" << rootSys.first[0] << ", " << rootSys.first[1] << ")" << std::endl;
    std::cout << "Iterations: " << rootSys.second << std::endl;
    std::cout << "FixedPoint solver test finished" << std::endl << std::endl;
}


void testRungeKutta()
{
    //u' = tu
    //[0, 1]
    //u(0) = 1
    //u = exp(0.5t^2)
    using Num::Ivp::RungeKuttaExplicit;
	using Num::Ivp::RungeKuttaImplicit;
    using Num::Ivp::ButcherTableau;
    using Num::Ivp::Methods;
    using Num::Ivp::Node;

	using Num::Equ::NeutonSystem;
	using Num::Ivp::DiffJacobian;

	using Num::Arg::MatNxM;
	using Num::Arg::VecN;



    auto solver = RungeKuttaExplicit<2, double, double>(Methods<double>::classic4());
    std::vector<RungeKuttaExplicit<2, double, double>::NextNode> nodesScalar;

    int n = 1000;
    double h = 1.0 / n;
    double arg = 0.0;
    double val = 1.0;
	
    for (int i = 0; i <= n; i++)
    {
        nodesScalar.push_back(
            solver.solve(
                [] (double t, double u) -> double
                {
                    return t * u;
                }
                , arg + h * i
                , val
                , h
            )
        );

        val = nodesScalar.back().second;
    }




    //D ^ 2(x)+2D(x) + x = 0
    //
    //x(t) = (At + B)exp(-t)
    //
    //[0, 1]
    //x(0) = +1
    //x'(0) = -1
    //x(t) = exp(-t)
    //
    // { x0' = x1,
    //{
    // { x1' = -2x1 - x0
    using Vec = VecN<double, 2>;
	using Mat = MatNxM<double, 2, 2>;


	auto func = [] (double t, Vec u) -> Vec
	{
		Vec res;

		const auto& x0 = u[0];
		const auto& x1 = u[1];

		res[0] = x1;
		res[1] = -2.0 * x1 - x0;

		return res;
	};

    auto solverSys4 = RungeKuttaExplicit<2, double, Vec>(Methods<double>::three_eighths_rule4());
    std::vector<decltype(solverSys4)::NextNode> nodesSystem4;


	auto tableau = Methods<double>::gaussLegendre6();
	auto solverSys6 = RungeKuttaImplicit<2, double, Vec, decltype(tableau)>(std::move(tableau), 1e-15, 20);
	std::vector<decltype(solverSys6)::NextNode> nodesSystem6;


    int nSys = 1000;
    double hSys = 1.0 / nSys;
    double argSys = 0.0;

    Vec    valSys4 = Vec(+1.0, -1.0);
	Vec    valSys6 = Vec(+1.0, -1.0);
    for (int i = 0; i <= nSys; i++)
    {
        nodesSystem4.push_back(
            solverSys4.solve(
                func
                , argSys + hSys * i
                , valSys4
                , hSys
            )
        );

		nodesSystem6.push_back(
			solverSys6.solve(
				func
				, DiffJacobian<Vec, Mat, decltype(func)>(std::move(func), 1e-8)
				, argSys + hSys * i
				, valSys6
				, hSys
			)
		);

        valSys4 = nodesSystem4.back().second;
		valSys6 = nodesSystem6.back().second;
    }


	



    std::cout.precision(15);
    std::cout << "RungeKuttaScalar test" << std::endl;

    std::cout << "Comparing numeric solution with exact" << std::endl;

    const auto&[t, y] = nodesScalar.back();

    std::cout << "Arg : " << t << std::endl;
    std::cout << "Exact  : " << exp(0.5 * t * t) << std::endl;
    std::cout << "Numeric: " << y << std::endl << std::endl;

    std::cout << "RungeKuttaScalar test finished" << std::endl << std::endl;


    std::cout << "RungeKuttaSystem test" << std::endl;

    std::cout << "Comparing numeric solution with exact" << std::endl;

	double tSys; Vec ySys;
    std::tie(tSys, ySys) = nodesSystem4.back();
	std::cout << "RungeKuttaExplicit:" << std::endl;
    std::cout << "Arg : " << tSys << std::endl;
    std::cout << "Exact  : " << "( " << exp(-tSys) << ", " << -exp(-tSys) << ")" << std::endl;
    std::cout << "Numeric: " << "( " << ySys[0]    << ", " << ySys[1]     << ")" << std::endl << std::endl;

	std::tie(tSys, ySys) = nodesSystem6.back();
	std::cout << "RungeKuttaImplicit:" << std::endl;
	std::cout << "Arg : " << tSys << std::endl;
	std::cout << "Exact  : " << "( " << exp(-tSys) << ", " << -exp(-tSys) << ")" << std::endl;
	std::cout << "Numeric: " << "( " << ySys[0]    << ", " << ySys[1]     << ")" << std::endl << std::endl;

    std::cout << "RungeKuttaSystem test finished" << std::endl << std::endl;
}


template<int N>
class GaussZeidel
{
public:
    using Mat = Num::Arg::MatNxM<double, N, N>;
    using Vec = Num::Arg::VecN<double, N>;

    GaussZeidel(Mat mat, Vec vec) : m_mat(mat), m_vec(vec)
    {}

    Vec operator() (const Vec& x)
    {
        Vec next;

        for (int i = 0; i < N; i++)
        {
            double sum = m_vec[i];

            for (int j = 0; j < i; j++)
            {
                sum -= m_mat[i][j] * next[j];
            }
            for (int j = i + 1; j < N; j++)
            {
                sum -= m_mat[i][j] * x[j];
            }

            next[i] = sum / m_mat[i][i];
        }

        return next;
    }


private:
    Mat m_mat;
    Vec m_vec;
};


int main()
{
    testLinAlg();
    testFixed();
    testNeuton();
    testRungeKutta();

	

    std::cin.get();
    return 0;
}