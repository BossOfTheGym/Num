#pragma once

#include "../LinAlg/Gauss.h"
#include "../Arg/ArgN.h"
#include "../Common/Common.h"
#include "../Utils/Utils.h"

#include "Common.h"



namespace Num
{
    namespace Equ
    {
        template<class Scalar, int N>
        using neu_sys_norm = Utils::Norm<Arg::VecN<Scalar, N>, Utils::Octo>;
        
        template<
              class ScalarType
            , int N
            , class Norm
            , template<class Scalar = ScalarType, int ROWS = N, int COLS = N> class MatrixType = Arg::MatNxM
            , template<class Scalar = ScalarType, int SIZE = N> class VectorType = Arg::VecN
            , class Solver = LinAlg::GaussEliminationSingle<ScalarType, N>
        >
        class NeutonSystem : public IterativeSolverBase<ScalarType>
        {
        public:
            using Base = IterativeSolverBase<ScalarType>;
            using Base::limit;
            using Base::eps;

            using Scalar = ScalarType;
            using Matrix = MatrixType<Scalar, N, N>;
            using Vector = VectorType<Scalar, N>;


		public:
            NeutonSystem() : Base(), m_solver(), m_norm()
            {}

            NeutonSystem(int iterationsLimit, Scalar eps = EPS<Scalar>) 
                : Base(iterationsLimit, eps)
                , m_norm()
            {}

            NeutonSystem(const NeutonSystem& ns) = default;
            NeutonSystem(NeutonSystem&& ns)      = default;


		public:
			template<class Function, class Jacobian>
            Vector solve(
                  Function&& func
                , Jacobian&& jacobian
                , const Vector& start
            )
            {
                Vector xPrev;
                Vector xNext;
                Vector delta;

                Matrix mat;

                int iterations = 0;


                xNext = start;
                do
                {
                    iterations++;

                    xPrev = xNext;

                    mat = jacobian(xPrev);
                    delta = -func(xPrev);

                    m_solver.solve(mat, delta);

                    xNext = xPrev + delta;
                } while (m_norm(delta) > eps() && iterations < limit());

                return xNext;
            }


        private:
            Solver m_solver;
            Norm   m_norm;
        };


        template<class Scalar>
        using neu_sc_norm = Utils::Norm<Arg::VecN<Scalar, 1>, Utils::Octo>;


        template<class Argument, class Norm>
        class NeutonScalar : public IterativeSolverBase<Argument>
        {
        public:
            using SystemSolver = NeutonSystem<Argument, 1, Norm>;

            using Base = IterativeSolverBase<Argument>;
            using Base::limit;
            using Base::eps;

            using Scalar = typename SystemSolver::Scalar;
            using Matrix = typename SystemSolver::Matrix;
            using Vector = typename SystemSolver::Vector;


        public:
            NeutonScalar() : Base(), m_solver()
            {}

            NeutonScalar(int iterationsLimit, Scalar eps = EPS<Scalar>) 
                : Base(iterationsLimit, eps)
                , m_solver(iterationsLimit, eps)
            {}

        public:
            template<class Function, class Derivative>
            Argument solve(Function&& func, Derivative&& deriv, const Argument& start)
            {
                auto wrappedFunc = [&] (const Vector& arg)
                {
                    return Vector(func(arg[0]));
                };

                auto wrappedJacob = [&] (const Vector& arg)
                {
                    return Matrix(deriv(arg[0]));
                };

                return m_solver.solve(wrappedFunc, wrappedJacob, Vector(start))[0];
            }

        private:
            SystemSolver m_solver;
        };
    }
}