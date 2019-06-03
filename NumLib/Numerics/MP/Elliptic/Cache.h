#ifndef Elliptic_Cache_h
#define Elliptic_Cache_h

#include "../Utility/Utility.h"



namespace Elliptic
{
	template<
		  class Cache

		, class Function

		, class LeftBoundFunction
		, class RightBoundFunction
		, class LowerBoundFunction
		, class UpperBoundFunction

		, class RightLowerCornerFunction
		, class RightUpperCornerFunction
		, class LeftLowerCornerFunction
		, class LeftUpperCornerFunction
	>
	void cacheEllipticFunctions( 
		int nt, int nx
		, Cache& cache

		, Function&&  function

		, LeftBoundFunction&&   leftBoundFunction
		, RightBoundFunction&& rightBoundFunction
		, LowerBoundFunction&& lowerBoundFunction
		, UpperBoundFunction&& upperBoundFunction

		, RightLowerCornerFunction&& rightLowerFunction
		, RightUpperCornerFunction&& rightUpperFunction
		, LeftLowerCornerFunction&&   leftLowerFunction
		, LeftUpperCornerFunction&&   leftUpperFunction 
	)
	{
		//lower
		cache[0][0] = leftLowerFunction();
		for (int j = 1; j < nx; j++)
		{
			cache[0][j] = lowerBoundFunction(j);
		}
		cache[0][nx] = rightLowerFunction();

		//inner
		for (int i = 1; i < nt; i++)
		{
			cache[i][0] = leftBoundFunction(i);
			for (int j = 1; j < nx; j++)
			{
				cache[i][j] = function(i, j);
			}
			cache[i][nx] = rightBoundFunction(i);
		}

		//upper
		cache[nt][0] = leftUpperFunction();
		for(int j = 1; j < nx; j++)
		{
			cache[nt][j] = upperBoundFunction(j);
		}
		cache[nt][nx] = rightUpperFunction();
	}



	template<
		  class Cache, class Solver, class Scalar

		, class Equation, class Function

		, class LeftBound , class LeftBoundFunction
		, class RightBound, class RightBoundFunction
		, class LowerBound, class LowerBoundFunction
		, class UpperBound, class UpperBoundFunction

		, class RightLowerCorner, class RightLowerFunction
		, class RightUpperCorner, class RightUpperFunction
		, class LeftLowerCorner , class LeftLowerFunction
		, class LeftUpperCorner , class LeftUpperFunction
	>
	auto solveCached(
		Cache& cache, Solver& solver

		, int nt, int nx, Scalar eps, int iterationLimit

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
		cacheEllipticFunctions(
			nt, nx, cache
			, function
			, leftBoundFunction
			, rightBoundFunction
			, lowerBoundFunction
			, upperBoundFunction
			, rightLowerFunction
			, rightUpperFunction
			, leftLowerFunction
			, leftUpperFunction
		);

		auto cachedFunction = [&] (int i, int j)
		{
			return cache[i][j];
		};


		auto cachedLeftBoundFunction = [&] (int i)
		{
			return cache[i][0];
		};

		auto cachedRightBoundFunction = [&] (int i)
		{
			return cache[i][nx];
		};

		auto cachedLowerBoundFunction = [&] (int j)
		{
			return cache[0][j];
		};

		auto cachedUpperBoundFunction = [&] (int j)
		{
			return cache[nt][j];
		};


		auto cachedRightLowerFunction = [&] ()
		{
			return cache[0][nx];
		};

		auto cachedRightUpperFunction = [&] ()
		{
			return cache[nt][nx];
		};

		auto cachedLeftLowerFunction = [&] ()
		{
			return cache[0][0];
		};

		auto cachedLeftUpperFunction = [&] ()
		{
			return cache[nt][0];
		};


		return solver.solve(
			nt, nx, eps, iterationLimit

			, equation        , cachedFunction

			, leftBound       , cachedLeftBoundFunction
			, rightBound      , cachedRightBoundFunction
			, lowerBound      , cachedLowerBoundFunction
			, upperBound      , cachedUpperBoundFunction

			, rightLowerCorner, cachedRightLowerFunction
			, rightUpperCorner, cachedRightUpperFunction
			, leftLowerCorner , cachedLeftLowerFunction
			, leftUpperCorner , cachedLeftUpperFunction
		);
	}
}

#endif