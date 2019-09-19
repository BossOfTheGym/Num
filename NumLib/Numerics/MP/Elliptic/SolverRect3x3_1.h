#ifndef SolverRect3x3_h
#define SolverRect3x3_h


#include "../Utility/Utility.h"

#include <utility>
#include <algorithm>
#include <memory>


namespace Elliptic
{	
	template<class Scalar>
	using Solution = Utility::Rect<Scalar>;


	
	//Solvers for rectangular region for template 3x3, 1 dimensional variable
	//
	// equation = function (all coefs are from the left, all values are from the right)
	// x2(j)
	// ^
	// |
	// |  (i + 1, j - 1)  (i + 1, j)  (i + 1, j + 1)        _____
	// |                                                  2|_____| next layer               
	// |  (i    , j - 1)  (i    , j)  (i    , j + 1) ---> 1|_____| current layer
	// |                                                  0|_____| previous layer
	// |  (i - 1, j - 1)  (i - 1, j)  (i + 1, j + 1)        0 1 2
	// |
	//  ---------------------------------------------> x1(i)
	//
	//class Equation        , class Function            - takes (i, j) -> Template / Value
	//
	//class LeftBound       , class LeftBoundFunction   - takes ( , j) -> Template / Value
	//class RightBound      , class RightBoundFunction  - takes ( , j) -> Template / Value
	//class LowerBound      , class LowerBoundFunction  - takes (i,  ) -> Template / Value
	//class UpperBound      , class UpperBoundFunction  - takes (i,  ) -> Template / Value
	//
	//class RightLowerCorner, class RightLowerFunction  - takes (void) -> Template / Value
	//class RightUpperCorner, class RightUpperFunction  - takes (void) -> Template / Value
	//class LeftLowerCorner , class LeftLowerFunction   - takes (void) -> Template / Value
	//class LeftUpperCorner , class LeftUpperFunction   - takes (void) -> Template / Value

	
	template<class Scalar>
	class SolverRect3x3_1_Simple
	{
	public:
		using SolutionT = Solution<Scalar>;


	public:
		template<
			class Equation        , class Function

			, class LeftBound       , class LeftBoundFunction
			, class RightBound      , class RightBoundFunction
			, class LowerBound      , class LowerBoundFunction
			, class UpperBound      , class UpperBoundFunction

			, class RightLowerCorner, class RightLowerFunction
			, class RightUpperCorner, class RightUpperFunction
			, class LeftLowerCorner , class LeftLowerFunction
			, class LeftUpperCorner , class LeftUpperFunction
		>
		SolutionT solve(
			int nx1, int nx2, Scalar eps, int iterationLimit
		
			, Equation&&                 equation, Function&&                     function
		
			, LeftBound&&               leftBound, LeftBoundFunction&&   leftBoundFunction
			, RightBound&&             rightBound, RightBoundFunction&& rightBoundFunction
			, LowerBound&&             lowerBound, LowerBoundFunction&& lowerBoundFunction
			, UpperBound&&             upperBound, UpperBoundFunction&& upperBoundFunction
		
			, RightLowerCorner&& rightLowerCorner, RightLowerFunction&& rightLowerFunction
			, RightUpperCorner&& rightUpperCorner, RightUpperFunction&& rightUpperFunction
			, LeftLowerCorner&&   leftLowerCorner, LeftLowerFunction&&   leftLowerFunction
			, LeftUpperCorner&&   leftUpperCorner, LeftUpperFunction&&   leftUpperFunction
		)
		{
			//common parameters
			auto prevSolution = std::make_unique<SolutionT>(nx1 + 1, nx2 + 1, 0.0);
			auto nextSolution = std::make_unique<SolutionT>(nx1 + 1, nx2 + 1, 0.0);

			int i = 0;

			int iterations = 0;


			// first layer : left lower corner + lower bound + right lower corner
			auto firstLayerStep = [&] ()
			{
				auto& prev = *prevSolution;
				auto& next = *nextSolution;

				//left lower corner
				auto llCorner      = leftLowerCorner();
				auto llCornerValue = leftLowerFunction();

				llCornerValue -=
					+ prev[1][0] * llCorner[2][1] + prev[1][1] * llCorner[2][2] 
					+                             + prev[0][1] * llCorner[1][2];

				next[0][0] = llCornerValue / llCorner[1][1];

				//lower bound
				for (int j = 1; j < nx2; j++)
				{
					auto lBound      = lowerBound(j);
					auto lBoundValue = lowerBoundFunction(j);

					lBoundValue -=
						+ prev[1][j - 1] * lBound[2][0] + prev[1][j] * lBound[2][1] + prev[1][j + 1] * lBound[2][2] 
						+ prev[0][j - 1] * lBound[1][0] +                           + prev[0][j + 1] * lBound[1][2];

					next[0][j] = lBoundValue / lBound[1][1];
				}

				//right lower corner
				auto rlCorner      = rightLowerCorner();
				auto rlCornerValue = rightLowerFunction();

				rlCornerValue -=
					+ prev[1][nx2 - 1] * rlCorner[2][0] + prev[1][nx2] * rlCorner[2][1] 
					+ prev[0][nx2 - 1] * rlCorner[1][0]                               ;

				next[0][nx2] = rlCornerValue / rlCorner[1][1];	
			};


			// inner layer : left bound + inner + right bound
			auto innerLayerStep = [&] ()
			{
				auto& prev = *prevSolution;
				auto& next = *nextSolution;

				//left bound
				auto lBound      = leftBound(i);
				auto lBoundValue = leftBoundFunction(i);

				lBoundValue -= 
					+ prev[i + 1][0] * lBound[2][1] + prev[i + 1][1] * lBound[2][2]
					+                               + prev[i    ][1] * lBound[1][2]
					+ prev[i - 1][0] * lBound[0][1] + prev[i - 1][1] * lBound[0][2];

				next[i][0] = lBoundValue / lBound[1][1];

				//inner
				for (int j = 1; j < nx2; j++)
				{
					auto inner      = equation(i, j);
					auto innerValue = function(i, j);

					innerValue -=
						+ prev[i + 1][j - 1] * inner[2][0] + prev[i + 1][j] * inner[2][1] + prev[i + 1][j + 1] * inner[2][2] 
						+ prev[i    ][j - 1] * inner[1][0] +                              + prev[i    ][j + 1] * inner[1][2] 
						+ prev[i - 1][j - 1] * inner[0][0] + prev[i - 1][j] * inner[0][1] + prev[i - 1][j + 1] * inner[0][2];

					next[i][j] = innerValue / inner[1][1];
				}
				 
				//right bound
				auto rBound      = rightBound(i);
				auto rBoundValue = rightBoundFunction(i);

				rBoundValue -=
					+ prev[i + 1][nx2 - 1] * rBound[2][0] + prev[i + 1][nx2] * rBound[2][1]
					+ prev[i    ][nx2 - 1] * rBound[1][0] +
					+ prev[i - 1][nx2 - 1] * rBound[0][0] + prev[i - 1][nx2] * rBound[0][1];

				next[i][nx2] = rBoundValue / rBound[1][1];
			};


			// last layer : left upper corner + upperBound + right upper corner
			auto lastLayerStep = [&] ()
			{
				auto& prev = *prevSolution;
				auto& next = *nextSolution;

				// left upper corner
				auto luCorner      = leftUpperCorner();
				auto luCornerValue = leftUpperFunction();

				luCornerValue -=
					+                                  + prev[nx1    ][1] * luCorner[1][2] 
					+ prev[nx1 - 1][0] * luCorner[0][1] + prev[nx1 - 1][1] * luCorner[0][2];

				next[nx1][0] = luCornerValue / luCorner[1][1];;

				// upper bound
				for (int j = 1; j < nx2; j++)
				{
					auto uBound      = upperBound(j);
					auto uBoundValue = upperBoundFunction(j);

					uBoundValue -=
						+ prev[nx1    ][j - 1] * uBound[1][0] +                                + prev[nx1    ][j + 1] * uBound[1][2] 
						+ prev[nx1 - 1][j - 1] * uBound[0][0] + prev[nx1 - 1][j] * uBound[0][1] + prev[nx1 - 1][j + 1] * uBound[0][2];

					next[nx1][j] = uBoundValue / uBound[1][1];
				}


				// right upper corner
				auto ruCorner      = rightUpperCorner();
				auto ruCornerValue = rightUpperFunction();

				ruCornerValue -=
					+ prev[nx1    ][nx2 - 1] * ruCorner[1][0] + 
					+ prev[nx1 - 1][nx2 - 1] * ruCorner[0][0] + prev[nx1 - 1][nx2] * ruCorner[0][1];

				next[nx1][nx2] = ruCornerValue  / ruCorner[1][1];;
			};

			//norm
			auto norm = [&] ()
			{
				auto& curr = *prevSolution;
				auto& next = *nextSolution;

				Scalar diff = static_cast<Scalar>(0);

				auto iter1 = curr.begin(); auto end = curr.end();
				auto iter2 = next.begin();
				while (iter1 != end)
				{
					diff = std::max(diff, std::abs(*iter1 - *iter2));

					++iter1;
					++iter2;
				}

				return diff;
			};


			//main loop
			do
			{
				++iterations;

				std::swap(prevSolution, nextSolution);

				firstLayerStep();
				for (i = 1; i < nx1; i++)
				{		
					innerLayerStep();
				}
				lastLayerStep();
			}
			while (norm() > eps && iterations < iterationLimit);

			return SolutionT(std::move(*nextSolution));
		}
	};



	template<class Scalar>
	class SolverRect3x3_1_Zeidel
	{
	public:
		using SolutionT = Solution<Scalar>;


	public:
		template<
			class Equation        , class Function

			, class LeftBound       , class LeftBoundFunction
			, class RightBound      , class RightBoundFunction
			, class LowerBound      , class LowerBoundFunction
			, class UpperBound      , class UpperBoundFunction

			, class RightLowerCorner, class RightLowerFunction
			, class RightUpperCorner, class RightUpperFunction
			, class LeftLowerCorner , class LeftLowerFunction
			, class LeftUpperCorner , class LeftUpperFunction
		>
		SolutionT solve(
			int nx1, int nx2, Scalar eps, int iterationLimit

			, Equation&&                 equation, Function&&                     function

			, LeftBound&&               leftBound, LeftBoundFunction&&   leftBoundFunction
			, RightBound&&             rightBound, RightBoundFunction&& rightBoundFunction
			, LowerBound&&             lowerBound, LowerBoundFunction&& lowerBoundFunction
			, UpperBound&&             upperBound, UpperBoundFunction&& upperBoundFunction

			, RightLowerCorner&& rightLowerCorner, RightLowerFunction&& rightLowerFunction
			, RightUpperCorner&& rightUpperCorner, RightUpperFunction&& rightUpperFunction
			, LeftLowerCorner&&   leftLowerCorner, LeftLowerFunction&&   leftLowerFunction
			, LeftUpperCorner&&   leftUpperCorner, LeftUpperFunction&&   leftUpperFunction
		)
		{
			//common parameters
			auto solution = SolutionT(nx1 + 1, nx2 + 1, 0.0);

			int i = 0;

			int iterations = 0;

			Scalar prev = static_cast<Scalar>(0);;
			Scalar norm = static_cast<Scalar>(0);


			// first layer : left lower corner + lower bound + right lower corner
			auto firstLayerStep = [&] ()
			{
				//left lower corner
				auto llCorner      = leftLowerCorner();
				auto llCornerValue = leftLowerFunction();

				prev = solution[0][0];

				llCornerValue -= 
					+ solution[1][0] * llCorner[2][1] + solution[1][1] * llCorner[2][2] 
					+                                 + solution[0][1] * llCorner[1][2];

				solution[0][0] = llCornerValue / llCorner[1][1];

				norm = std::max(norm, std::abs(prev - solution[0][0]));


				//lower bound
				for (int j = 1; j < nx2; j++)
				{
					auto lBound      = lowerBound(j);
					auto lBoundValue = lowerBoundFunction(j);

					prev = solution[0][j];

					lBoundValue -= 
						+ solution[1][j - 1] * lBound[2][0] + solution[1][j] * lBound[2][1] + solution[1][j + 1] * lBound[2][2]
						+ solution[0][j - 1] * lBound[1][0] +                               + solution[0][j + 1] * lBound[1][2];

					solution[0][j] = lBoundValue / lBound[1][1];

					norm = std::max(norm, std::abs(prev - solution[0][j]));
				}


				//right lower corner
				auto rlCorner      = rightLowerCorner();
				auto rlCornerValue = rightLowerFunction();

				prev = solution[0][nx2];

				rlCornerValue -=
					+ solution[1][nx2 - 1] * rlCorner[2][0] + solution[1][nx2] * rlCorner[2][1] 
					+ solution[0][nx2 - 1] * rlCorner[1][0];

				solution[0][nx2] = rlCornerValue  / rlCorner[1][1];

				norm = std::max(norm, std::abs(prev - solution[0][nx2]));
			};


			// inner layer : left bound + inner + right bound
			auto innerLayerStep = [&] ()
			{
				//left bound
				auto lBound      = leftBound(i);
				auto lBoundValue = leftBoundFunction(i);

				prev = solution[i][0];

				lBoundValue -= 
					+ solution[i + 1][0] * lBound[2][1] + solution[i + 1][1] * lBound[2][2]
					+                                   + solution[i    ][1] * lBound[1][2]
					+ solution[i - 1][0] * lBound[0][1] + solution[i - 1][1] * lBound[0][2];

				solution[i][0] = lBoundValue / lBound[1][1];

				norm = std::max(norm, std::abs(prev - solution[i][0]));


				//inner
				for (int j = 1; j < nx2; j++)
				{
					auto inner      = equation(i, j);
					auto innerValue = function(i, j);

					prev = solution[i][j];

					innerValue -=
						+ solution[i + 1][j - 1] * inner[2][0] + solution[i + 1][j] * inner[2][1] + solution[i + 1][j + 1] * inner[2][2] 
						+ solution[i    ][j - 1] * inner[1][0] +                                  + solution[i    ][j + 1] * inner[1][2] 
						+ solution[i - 1][j - 1] * inner[0][0] + solution[i - 1][j] * inner[0][1] + solution[i - 1][j + 1] * inner[0][2];

					solution[i][j] = innerValue / inner[1][1];

					norm = std::max(norm, std::abs(prev - solution[i][j]));
				}


				//right bound
				auto rBound      = rightBound(i);
				auto rBoundValue = rightBoundFunction(i);

				prev = solution[i][nx2];

				rBoundValue -= 
					+ solution[i + 1][nx2 - 1] * rBound[2][0] + solution[i + 1][nx2] * rBound[2][1]
					+ solution[i    ][nx2 - 1] * rBound[1][0] + 
					+ solution[i - 1][nx2 - 1] * rBound[0][0] + solution[i - 1][nx2] * rBound[0][1];

				solution[i][nx2] = rBoundValue / rBound[1][1];

				norm = std::max(norm, std::abs(prev - solution[i][nx2]));
			};


			// last layer : left upper corner + upperBound + right upper corner
			auto lastLayerStep = [&] ()
			{
				// left upper corner
				auto luCorner      = leftUpperCorner();
				auto luCornerValue = leftUpperFunction();

				prev = solution[nx1][0];

				luCornerValue -=
					+ solution[nx1    ][1] * luCorner[1][2] 
					+ solution[nx1 - 1][0] * luCorner[0][1] + solution[nx1 - 1][1] * luCorner[0][2];

				solution[nx1][0] = luCornerValue / luCorner[1][1];

				norm = std::max(norm, std::abs(prev - solution[nx1][0]));


				// upper bound
				for (int j = 1; j < nx2; j++)
				{
					auto uBound      = upperBound(j);
					auto uBoundValue = upperBoundFunction(j);

					prev = solution[nx1][j];

					uBoundValue -= 
						+ solution[nx1    ][j - 1] * uBound[1][0] +                                    + solution[nx1    ][j + 1] * uBound[1][2] 
						+ solution[nx1 - 1][j - 1] * uBound[0][0] + solution[nx1 - 1][j] * uBound[0][1] + solution[nx1 - 1][j + 1] * uBound[0][2];

					solution[nx1][j] = uBoundValue / uBound[1][1];

					norm = std::max(norm, std::abs(prev - solution[nx1][j]));
				}


				// right upper corner
				auto ruCorner      = rightUpperCorner();
				auto ruCornerValue = rightUpperFunction();

				prev = solution[nx1][nx2];

				ruCornerValue -= 
					+ solution[nx1    ][nx2 - 1] * ruCorner[1][0] + 
					+ solution[nx1 - 1][nx2 - 1] * ruCorner[0][0] + solution[nx1 - 1][nx2] * ruCorner[0][1];

				solution[nx1][nx2] = ruCornerValue / ruCorner[1][1];

				norm = std::max(norm, std::abs(prev - solution[nx1][nx2]));
			};


			//main loop
			do
			{
				++iterations;

				norm = static_cast<Scalar>(0);

				firstLayerStep();
				for (i = 1; i < nx1; i++)
				{		
					innerLayerStep();
				}
				lastLayerStep();
			}
			while (norm > eps && iterations < iterationLimit);

			//TODO: remove this
			std::cout << iterations << std::endl;

			return std::move(solution);
		}
	};



	template<class Scalar>
	class SolverRect3x3_1_Relaxation
	{
	public:
		using SolutionT = Solution<Scalar>;


	public:
		SolverRect3x3_1_Relaxation(Scalar relaxationParameter) : m_relax(relaxationParameter)
		{}

		SolverRect3x3_1_Relaxation(const SolverRect3x3_1_Relaxation& sr) = default;
		SolverRect3x3_1_Relaxation(SolverRect3x3_1_Relaxation&& sr)      = default;

		~SolverRect3x3_1_Relaxation() = default;

		SolverRect3x3_1_Relaxation& operator = (const SolverRect3x3_1_Relaxation& sr) = default;
		SolverRect3x3_1_Relaxation& operator = (SolverRect3x3_1_Relaxation&& sr)      = default;


	public:
		template<
			class Equation        , class Function

			, class LeftBound       , class LeftBoundFunction
			, class RightBound      , class RightBoundFunction
			, class LowerBound      , class LowerBoundFunction
			, class UpperBound      , class UpperBoundFunction

			, class RightLowerCorner, class RightLowerFunction
			, class RightUpperCorner, class RightUpperFunction
			, class LeftLowerCorner , class LeftLowerFunction
			, class LeftUpperCorner , class LeftUpperFunction
		>
		SolutionT solve(
			int nx1, int nx2, Scalar eps, int iterationLimit

			, Equation&&                 equation, Function&&                     function

			, LeftBound&&               leftBound, LeftBoundFunction&&   leftBoundFunction
			, RightBound&&             rightBound, RightBoundFunction&& rightBoundFunction
			, LowerBound&&             lowerBound, LowerBoundFunction&& lowerBoundFunction
			, UpperBound&&             upperBound, UpperBoundFunction&& upperBoundFunction

			, RightLowerCorner&& rightLowerCorner, RightLowerFunction&& rightLowerFunction
			, RightUpperCorner&& rightUpperCorner, RightUpperFunction&& rightUpperFunction
			, LeftLowerCorner&&   leftLowerCorner, LeftLowerFunction&&   leftLowerFunction
			, LeftUpperCorner&&   leftUpperCorner, LeftUpperFunction&&   leftUpperFunction
		)
		{
			//common parameters
			auto solution = SolutionT(nx1 + 1, nx2 + 1, 0.0);

			int i = 0;

			int iterations = 0;

			Scalar prev = static_cast<Scalar>(0);
			Scalar norm = static_cast<Scalar>(0);

			const Scalar w = m_relax;


			// first layer : left lower corner + lower bound + right lower corner
			auto firstLayerStep = [&] ()
			{
				//left lower corner
				auto llCorner      = leftLowerCorner();
				auto llCornerValue = leftLowerFunction();

				prev = solution[0][0];

				llCornerValue -= 
					+ solution[1][0] * llCorner[2][1] + solution[1][1] * llCorner[2][2] 
					+                             + solution[0][1] * llCorner[1][2];

				solution[0][0] = llCornerValue * w / llCorner[1][1] + solution[0][0] * (static_cast<Scalar>(1) - w);

				norm = std::max(norm, std::abs(prev - solution[0][0]));


				//lower bound
				for (int j = 1; j < nx2; j++)
				{
					auto lBound      = lowerBound(j);
					auto lBoundValue = lowerBoundFunction(j);

					prev = solution[0][j];

					lBoundValue -= 
						+ solution[1][j - 1] * lBound[2][0] + solution[1][j] * lBound[2][1] + solution[1][j + 1] * lBound[2][2]
						+ solution[0][j - 1] * lBound[1][0] +                           + solution[0][j + 1] * lBound[1][2];

					solution[0][j] = lBoundValue * w / lBound[1][1] + solution[0][j] * (static_cast<Scalar>(1) - w);

					norm = std::max(norm, std::abs(prev - solution[0][j]));
				}


				//right lower corner
				auto rlCorner      = rightLowerCorner();
				auto rlCornerValue = rightLowerFunction();

				prev = solution[0][nx2];

				rlCornerValue -=
					+ solution[1][nx2 - 1] * rlCorner[2][0] + solution[1][nx2] * rlCorner[2][1] 
					+ solution[0][nx2 - 1] * rlCorner[1][0];

				solution[0][nx2] = rlCornerValue * w / rlCorner[1][1] + solution[0][nx2] * (static_cast<Scalar>(1) - w);

				norm = std::max(norm, std::abs(prev - solution[0][nx2]));
			};


			// inner layer : left bound + inner + right bound
			auto innerLayerStep = [&] ()
			{
				//left bound
				auto lBound      = leftBound(i);
				auto lBoundValue = leftBoundFunction(i);

				prev = solution[i][0];

				lBoundValue -= 
					+ solution[i + 1][0] * lBound[2][1] + solution[i + 1][1] * lBound[2][2]
					+                               + solution[i    ][1] * lBound[1][2]
					+ solution[i - 1][0] * lBound[0][1] + solution[i - 1][1] * lBound[0][2];

				solution[i][0] = lBoundValue * w / lBound[1][1] + solution[i][0] * (static_cast<Scalar>(1) - w);

				norm = std::max(norm, std::abs(prev - solution[i][0]));


				//inner
				for (int j = 1; j < nx2; j++)
				{
					auto inner      = equation(i, j);
					auto innerValue = function(i, j);

					prev = solution[i][j];

					innerValue -=
						+ solution[i + 1][j - 1] * inner[2][0] + solution[i + 1][j] * inner[2][1] + solution[i + 1][j + 1] * inner[2][2] 
						+ solution[i    ][j - 1] * inner[1][0] +                                  + solution[i    ][j + 1] * inner[1][2] 
						+ solution[i - 1][j - 1] * inner[0][0] + solution[i - 1][j] * inner[0][1] + solution[i - 1][j + 1] * inner[0][2];

					solution[i][j] = innerValue * w / inner[1][1] + solution[i][j] * (static_cast<Scalar>(1) - w);

					norm = std::max(norm, std::abs(prev - solution[i][j]));
				}


				//right bound
				auto rBound      = rightBound(i);
				auto rBoundValue = rightBoundFunction(i);

				prev = solution[i][nx2];

				rBoundValue -= 
					+ solution[i + 1][nx2 - 1] * rBound[2][0] + solution[i + 1][nx2] * rBound[2][1]
					+ solution[i    ][nx2 - 1] * rBound[1][0] + 
					+ solution[i - 1][nx2 - 1] * rBound[0][0] + solution[i - 1][nx2] * rBound[0][1];

				solution[i][nx2] = rBoundValue * w / rBound[1][1] + solution[i][nx2] * (static_cast<Scalar>(1) - w);

				norm = std::max(norm, std::abs(prev - solution[i][nx2]));
			};


			// last layer : left upper corner + upperBound + right upper corner
			auto lastLayerStep = [&] ()
			{
				// left upper corner
				auto luCorner      = leftUpperCorner();
				auto luCornerValue = leftUpperFunction();

				prev = solution[nx1][0];

				luCornerValue -=
					+ solution[nx1    ][1] * luCorner[1][2] 
					+ solution[nx1 - 1][0] * luCorner[0][1] + solution[nx1 - 1][1] * luCorner[0][2];

				solution[nx1][0] = luCornerValue * w / luCorner[1][1] + solution[nx1][0] * (static_cast<Scalar>(1) - w);

				norm = std::max(norm, std::abs(prev - solution[nx1][0]));

				// upper bound
				for (int j = 1; j < nx2; j++)
				{
					auto uBound      = upperBound(j);
					auto uBoundValue = upperBoundFunction(j);

					prev = solution[nx1][j];

					uBoundValue -= 
						+ solution[nx1    ][j - 1] * uBound[1][0] +                                    + solution[nx1    ][j + 1] * uBound[1][2] 
						+ solution[nx1 - 1][j - 1] * uBound[0][0] + solution[nx1 - 1][j] * uBound[0][1] + solution[nx1 - 1][j + 1] * uBound[0][2];

					solution[nx1][j] = uBoundValue * w / uBound[1][1] + solution[nx1][j] * (static_cast<Scalar>(1) - w);

					norm = std::max(norm, std::abs(prev - solution[nx1][j]));
				}


				// right upper corner
				auto ruCorner      = rightUpperCorner();
				auto ruCornerValue = rightUpperFunction();

				prev = solution[nx1][nx2];

				ruCornerValue -= 
					+ solution[nx1    ][nx2 - 1] * ruCorner[1][0] + 
					+ solution[nx1 - 1][nx2 - 1] * ruCorner[0][0] + solution[nx1 - 1][nx2] * ruCorner[0][1];

				solution[nx1][nx2] = ruCornerValue * w / ruCorner[1][1] + solution[nx1][nx2] * (static_cast<Scalar>(1) - w);

				norm = std::max(norm, std::abs(prev - solution[nx1][nx2]));
			};


			//main loop
			do
			{
				++iterations;

				//TODO: remove this
				if (iterations % 1000 == 0)
				{
					std::cout << iterations << std::endl;
				}


				norm = static_cast<Scalar>(0);

				firstLayerStep();
				for (i = 1; i < nx1; i++)
				{		
					innerLayerStep();
				}
				lastLayerStep();
			}
			while (norm > eps && iterations < iterationLimit);

			//TODO: and this
			std::cout << iterations << std::endl;

			return std::move(solution);
		}


	private:
		const Scalar m_relax;
	};
}

#endif