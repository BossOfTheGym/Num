#ifndef EQU_UTILS_H
#define EQU_UTILS_H

#include <utility>
#include <type_traits>


namespace Num
{
	namespace Equ
	{
		template<class Value, class Function>
		class DiffDerivative : public Function
		{
		public:
			template<class FuncType, class ValueType>
			DiffDerivative(FuncType&& function, ValueType&& eps) 
				: Function(std::forward<FuncType>(function))
				, mEps(std::forward<ValueType>(eps))
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

		template<class Function, class Value>
		auto make_diff_derivative(Function&& function, Value&& eps)
		{
			return DiffDerivative<
				  std::remove_cv_t<std::remove_reference_t<Value>>
				, std::remove_cv_t<std::remove_reference_t<Function>>
			>(std::forward<Function>(function), std::forward<Value>(eps));
		}



		template< class Vector, class Matrix, class Function>
		class DiffJacobian : public Function
		{
			static_assert(Vector::SIZE == Matrix::COLS && Vector::SIZE == Matrix::ROWS, "Non-matching sizes");

		public:
			using Scalar = typename Vector::Scalar;

		public:
			template<class FuncType, class ScalarType>
			DiffJacobian(FuncType&& function, ScalarType&& eps) 
				: Function(std::forward<FuncType>(function))
				, mEps(std::forward<ScalarType>(eps))
			{}

			DiffJacobian(const DiffJacobian&) = default;
			DiffJacobian(DiffJacobian&&)      = default;


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

		template<class Vector, class Matrix, class Function, class Eps>
		auto make_diff_jacobian(Function&& function, Eps&& eps)
		{
			return DiffJacobian<
				  Vector
				, Matrix
				, std::remove_reference_t<std::remove_cv_t<Function>>
			>(std::forward<Function>(function), std::forward<Eps>(eps));
		}
	}
}

#endif