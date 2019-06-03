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
        template<class Value>
        class NeutonScalar : public IterativeSolverBase<Value>
        {
        public:
            using Base = IterativeSolverBase<Value>;
            using Base::Base;
            using Base::limit;
            using Base::eps;


			template<class Function, class Derivative>
            Root<Value> solve(
                  Function&& func
                , Derivative&& deriv
                , const Value& start
            )
            {
				Value xPrev;
				Value xNext;

                int iterations = 0;

                xNext = start;
                do
                {
                    iterations++;

                    xPrev = xNext;

                    xNext = xPrev - func(xPrev) / deriv(xPrev);
                } while (abs(xNext - xPrev) > eps() && iterations < limit());

                return {xNext, iterations};
            }
        };


		template<class Value, class Function>
		class DiffDerivative : public Function
		{
		public:
			DiffDerivative(Function&& function, Value eps) 
				: Function(std::forward<Function>(function))
				, mEps(eps)
			{}

			DiffDerivative(const DiffDerivative& deriv) = default;
			DiffDerivative(DiffDerivative&& deriv)      = default;


		public:
			Value operator()(const Value& arg)
			{
				return (Function::operator()(arg + mEps) - Function::operator()(arg - mEps)) / (2 * mEps);
			}

		private:
			Value mEps;
		};
		


        template<
              class ScalarType
            , int N
            , template<class Scalar = ScalarType, int ROWS = N, int COLS = N> class MatrixType = Arg::MatNxM
            , template<class Scalar = ScalarType, int SIZE = N> class VectorType = Arg::VecN
            , class Solver = LinAlg::GaussEliminationSingle<ScalarType, N>
            , class Norm   = DefaultNorm<VectorType<ScalarType, N>>
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


            NeutonSystem() : Base(), m_solver(), m_norm()
            {}

            NeutonSystem(int iterationsLimit, Scalar eps = EPS<Scalar>) 
                : Base(iterationsLimit, eps)
                , m_norm()
            {}

            NeutonSystem(const NeutonSystem& ns) = default;
            NeutonSystem(NeutonSystem&& ns)      = default;



			template<class Function, class Jacobian>
            Root<Vector> solve(
                  Function&& func
                , Jacobian&& jacobian
                , Vector&& start
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

                return {xNext, iterations};
            }

			Norm& getNorm()
			{
				return m_norm;
			}


        private:
            Solver m_solver;
            Norm   m_norm;
        };


		template< 
			  class Vector
			, class Matrix
			, class Function
		>
		class DiffJacobian : public Function
		{
			static_assert(Vector::SIZE == Matrix::COLS && Vector::SIZE == Matrix::ROWS, "Non-matching sizes");

		public:
			using Scalar = typename Vector::Scalar;

		public:
			DiffJacobian(Function&& function, Scalar&& eps) 
				: Function(std::forward<Function>(function))
				, mEps(eps)
			{}

			DiffJacobian(const DiffJacobian& dj) = default;
			DiffJacobian(DiffJacobian&& dj)      = default;


			Matrix operator() (const Vector& arg)
			{
				const int N = Vector::SIZE;

				Matrix mat;

				for(int i = 0; i < N; i++)
				{
					Vector rightArg(arg);
					Vector  leftArg(arg);
					rightArg[i] += mEps;
					leftArg [i] -= mEps;
					

					Vector right = Function::operator()(rightArg);
					Vector left  = Function::operator()(leftArg);
					for (int j = 0; j < N; j++)
					{
						mat[j][i] = (right[j] - left[j]) / (2 * mEps);
					}
				}

				return mat;
			}


		private:
			Scalar mEps;
		};
    }
}