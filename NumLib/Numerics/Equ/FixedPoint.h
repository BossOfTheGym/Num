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
        template<
              class ScalarType
            , int N
            , class Norm
            , template<class Scalar = ScalarType, int ROWS = N, int COLS = N> class MatrixType = Arg::MatNxM
            , template<class Scalar = ScalarType, int SIZE = N> class VectorType = Arg::VecN
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
            Vector solve(
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

                return xNext;
            }

        private:
            Norm m_norm;
        };


        template<class Argument, class Norm>
        class FixedPointScalar : public IterativeSolverBase<Argument>
        {
        public:
            using SystemSolver = FixedPointSystem<Argument, 1, Norm>;

            using Base = IterativeSolverBase<Argument>;
            using Base::limit;
            using Base::eps;

            using Scalar = typename SystemSolver::ScalarType;
            using Matrix = typename SystemSolver::MatrixType;
            using Vector = typename SystemSolver::VectorType;


        public:
            FixedPointScalar() : Base(), m_solver()
            {}

            FixedPointScalar(int iterationsLimit, Scalar eps = EPS<Scalar>) 
                : Base(iterationsLimit, eps)
                , m_solver(iterationsLimit, eps)
            {}

        public:
            template<class Function>
            Argument solve(Function&& mapping, const Argument& start)
            {
                auto wrapped = [&] (const Vector& arg)
                {
                    return Vector(mapping(arg[0]));
                };

                return m_solver.solve(wrapped, Vector(start))[0];
            }

        private:
            SystemSolver m_solver;
        };
    }
}