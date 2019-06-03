#pragma once

#include <cmath>
#include <utility>

#include "../Common/Common.h"


namespace Num
{
    namespace Equ
    {
        template<class Value>
        using Root = std::pair<Value, int>;


        const int DEFAULT_ITERATIONS_LIMIT = 10;


        template<class Vector>
        using DefaultNorm = Utils::Norm<Vector, Utils::NormType::OCTO>;


        template<class Value>
        class IterativeSolverBase
        {
        public:
            IterativeSolverBase() : IterativeSolverBase(DEFAULT_ITERATIONS_LIMIT, EPS<Value>)
            {}

            IterativeSolverBase(int iterationsLimit, Value eps = EPS<Value>)
                : m_iterationsLimit(iterationsLimit)
                , m_eps(eps)
            {}

            IterativeSolverBase(const IterativeSolverBase& isb) = default;

            IterativeSolverBase(IterativeSolverBase&& isb) = default;



            int& limit()
            {
                return m_iterationsLimit;
            }

            const int& limit() const
            {
                return m_iterationsLimit;
            }

			Value& eps()
            {
                return m_eps;
            }

            const Value& eps() const
            {
                return m_eps;
            }


        private:
            int m_iterationsLimit;
			Value m_eps;
        };
    }
}