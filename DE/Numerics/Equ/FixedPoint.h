#pragma once

#include "../Arg/ArgN.h"
#include "../Common/Common.h"
#include "../Utils/Utils.h"

#include "Common.h"


//TODO : add different mappings-methods as functors
namespace Num
{
    namespace Equ
    {
        template<class Argument, class Value = Argument>
        class FixedPointScalar : public IterativeSolverBase<Value>
        {
        public:
            using Base = IterativeSolverBase<Value>;
            using Base::Base;
            using Base::limit;
            using Base::eps;


			template<class Function>
            Root<Argument> solve(
                  Function&& mapping
                , const Argument& start
            )
            {
                Argument xPrev;
                Argument xNext;
                
                int iterations = 0;

                xNext = start;
                do
                {
                    iterations++;

                    xPrev = xNext;

                    xNext = mapping(xPrev);
                } while (abs(xPrev - xNext) > eps() && iterations < limit());

                return {xNext, iterations};
            }
        };


        template<
              class ScalarType
            , int N
            , template<class Scalar = ScalarType, int ROWS = N, int COLS = N> class MatrixType = Arg::MatNxM
            , template<class Scalar = ScalarType, int SIZE = N> class VectorType = Arg::VecN
            , class Norm = DefaultNorm<VectorType<ScalarType, N>>
        >
        class FixedPointSystem : public IterativeSolverBase<ScalarType>
        {
        public:
            using Base = IterativeSolverBase<ScalarType>;
            using Base::limit;
            using Base::eps;

            using Scalar = ScalarType;
            using Matrix = MatrixType<Scalar, N, N>;
            using Vector = VectorType<Scalar, N>;


            FixedPointSystem() : Base(), m_norm()
            {}

            FixedPointSystem(int iterationsLimit, Scalar eps = EPS<Scalar>) 
                : Base(iterationsLimit, eps)
                , m_norm()
            {}

            FixedPointSystem(const FixedPointSystem& ns) = default;

            FixedPointSystem(FixedPointSystem&& ns) = default;

			template<class Function>
            Root<Vector> solve(
                  Function&& mapping
                , const Vector& start
            )
            {
                Vector xPrev;
                Vector xNext;

                int iterations = 0;

                xNext = start;
                do
                {
                    iterations++;

                    xPrev = xNext;

                    xNext = mapping(xPrev);
                } while (m_norm(xPrev - xNext) > eps() && iterations < limit());

                return {xNext, iterations};
            }

        private:
            Norm m_norm;
        };
    }
}