#ifndef EQU_UTILS_H
#define EQU_UTILS_H

#include <utility>
#include <type_traits>


namespace Num
{
	namespace Equ
	{
		template< class Vector, class Matrix, class Function>
		class DiffJacobian
		{
			static_assert(Vector::SIZE == Matrix::COLS && Vector::SIZE == Matrix::ROWS, "Non-matching sizes");

		public:
			using Scalar = typename Vector::Scalar;

		public:
			template<class FuncType, class ScalarType>
			DiffJacobian(FuncType&& function, ScalarType&& eps) 
				: m_func(std::forward<FuncType>(function))
				, m_eps(std::forward<ScalarType>(eps))
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
					rightArg[i] += m_eps;
					leftArg [i] -= m_eps;


					Vector right = m_func(rightArg);
					Vector left  = m_func(leftArg);
					for (int j = 0; j < N; j++)
					{
						mat[j][i] = (right[j] - left[j]) / (2 * m_eps);
					}
				}

				return mat;
			}


		private:
			Function m_func;
			Scalar m_eps;
		};

		template<class Vector, class Matrix, class Function, class Eps>
		auto make_diff_jacobian(Function&& function, Eps&& eps)
		{
			using FunctionType = std::remove_reference_t<std::remove_cv_t<Function>>;

			return DiffJacobian<Vector, Matrix, FunctionType>(
				std::forward<Function>(function), std::forward<Eps>(eps)
			);
		}

		template<class Function, class Scalar>
		auto make_diff_derivative(Function&& function, Scalar&& eps)
		{
			using Vec = Num::Arg::VecN<Scalar, 1>;
			using Mat = Num::Arg::MatNxM<Scalar, 1, 1>;

			auto wrapper = [&] (const Vec& arg)
			{
				return Vec(function(arg[0]));
			};
			
			return make_diff_jacobian<Vec, Mat>(wrapper, std::forward<Scalar>(eps));
		}
	}
}

#endif