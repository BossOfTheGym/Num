#ifndef SolverRect2x3_h
#define SolverRect2x3_h


#include "../Utility/Utility.h"

#include <utility>
#include <type_traits>

#include <array>
#include <vector>


namespace Parabolic
{
	//Solver for rectangular region for template 2x3
	// 1 dimensional variable
	//
	// t(i)
	// ^
	// |                                                     _____
	// |  (i    , j - 1)  (i    , j)  (i    , j + 1)       1|_____| current layer
	// |                                              ---> 0|_____| previous layer
	// |  (i - 1, j - 1)  (i - 1, j)  (i - 1, j + 1)         0 1 2
	// |
	//  ---------------------------------------------> x(j)
	template<class Scalar>
	class SolverRect2x3_1
	{
	public:
		using Solution = Utility::Rect<Scalar>;
		using Line = std::array<Scalar, 3>;


	public:
		template<
			  class Equation  , class Function
			, class LeftBound , class LeftBoundFunction
			, class RightBound, class RightBoundFunction
			,                   class LowerBoundFunction
		>
		Solution solve(
			  int nt, int nx
			, Equation&&   equation  , Function&&           function           // equation + value
			, LeftBound&&  leftBound , LeftBoundFunction&&  leftBoundFunction  // equation + value
			, RightBound&& rightBound, RightBoundFunction&& rightBoundFunction // equation + value
			,                          LowerBoundFunction&& lowerBoundFunction // value only
		)
		{
			using Template    = std::invoke_result_t<std::remove_reference_t<Equation>, int, int>;
			using SystemBlock = std::vector<Line>;


			// common parameters
			Solution solution(nt + 1, nx + 1, 0.0);
			SystemBlock block(nx + 1);



			//initialize lower bound
			auto initializeLowerBound = [&] () -> void
			{
				for (int j = 0; j <= nx; j++)
				{
					solution[0][j] = lowerBoundFunction(j);
				}
			};


			////initialize current block
			//reduce previous layer
			auto reducePreviousLayerFirst = [&] (int i, auto& coefs)
			{
				solution[i][0] -= solution[i - 1][0] * coefs[1];
				solution[i][0] -= solution[i - 1][1] * coefs[2];
			};

			auto reducePreviousLayerInner = [&] (int i, int j, auto& coefs)
			{
				solution[i][j] -= solution[i - 1][j - 1] * coefs[0];
				solution[i][j] -= solution[i - 1][j    ] * coefs[1];
				solution[i][j] -= solution[i - 1][j + 1] * coefs[2];
			};

			auto reducePreviousLayerLast = [&] (int i, auto& coefs)
			{
				solution[i][nx] -= solution[i - 1][nx - 1] * coefs[0];
				solution[i][nx] -= solution[i - 1][nx    ] * coefs[1];
			};

			
			auto initializeCurrentSystem = [&] (int i)
			{
				//left bound
				auto rb = rightBound(i);
				block[0][0] = rb[1][0];
				block[0][1] = rb[1][1];
				block[0][2] = rb[1][2];

				reducePreviousLayerFirst(i, rb[0]);

				//inner
				for (int j = 1; j < nx; j++)
				{
					auto inner = equation(i, j);
					block[j][0] = inner[1][0];
					block[j][1] = inner[1][1];
					block[j][2] = inner[1][2];

					reducePreviousLayerInner(i, j, inner[0]);
				}

				//right bound
				auto lb = leftBound(i);
				block[nx][0] = lb[1][0];
				block[nx][1] = lb[1][1];
				block[nx][2] = lb[1][2];

				reducePreviousLayerLast(i, lb[0]);
			};

			auto initializeCurrentConditions = [&] (int i)
			{
				solution[i][0] = leftBoundFunction(i);
				
				for (int j = 1; j < nx; j++)
				{
					solution[i][j] = function(i, j);
				}

				solution[i][nx] = rightBoundFunction(i);
			};

			auto initializeCurrentBlock = [&] (int i)
			{
				initializeCurrentConditions(i);
				initializeCurrentSystem(i);
			};


			//forward pass
			auto forwardPass = [&] (int i)
			{
				for (int j = 1; j <= nx; j++)
				{
					auto&[a0, b0, c0] = block[j - 1];
					auto&[a1, b1, c1] = block[j    ];

					solution[i][j - 1] /= b0;
					c0 /= b0;
					b0 /= b0;

					solution[i][j] -= a1 * solution[i][j - 1];
					b1 -= a1 * c0;
					a1 -= a1;
				}

				auto& b1 = block[nx][1];
				solution[i][nx] /= b1;
				b1 /= b1;
			};

			//backward pass
			auto backwardPass = [&] (int i)
			{
				for (int j = nx; j > 0; j--)
				{
					auto&[a, b, c] = block[j - 1];

					solution[i][j - 1] -= c * solution[i][j];
				}
			};


			//algorithm
			//1.
			initializeLowerBound();
			//2.
			for (int i = 1; i <= nt; i++)
			{
				initializeCurrentBlock(i);

				forwardPass(i);

				backwardPass(i);
			}
			//3.
			return solution;
		}
	};
}

#endif
