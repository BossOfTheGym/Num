#pragma once

#include <utility>
#include <type_traits>


#include "../Equ/Neuton.h"
#include "../Equ/Utils.h"

namespace Num
{
	namespace Ivp
	{
		template< class Vector, class Matrix, class Function>
		class DiffJacobian
		{
		public:
			static_assert(Vector::SIZE == Matrix::COLS && Vector::SIZE == Matrix::ROWS, "Non-matching sizes");

			using Scalar = typename Vector::Scalar;


		public:
			template<class FuncType, class EpsType>
			DiffJacobian(FuncType&& function, EpsType&& eps) 
				: m_func(std::forward<FuncType>(function))
				, m_eps(std::forward<EpsType>(eps))
			{}

			DiffJacobian(const DiffJacobian&) = default;
			DiffJacobian(DiffJacobian&&)      = default;


		public:
			Matrix operator() (const Scalar& arg, const Vector& val)
			{
				const int N = Vector::SIZE;

				auto bound = [&] (const Vector& value)
				{
					return m_func(arg, value);
				};

				return Equ::make_diff_jacobian<Vector, Matrix>(bound, m_eps)(val);
			}

		private:
			Function m_func;
			Scalar   m_eps;
		};


		template<class Vector, class Matrix, class Eps, class Function>
		auto make_diff_jacobian(Function&& function, Eps&& eps)
		{
			using FunctionType = std::remove_reference_t<std::remove_cv_t<Function>>;

			return DiffJacobian<Vector, Matrix, FunctionType>(
				std::forward<Function>(function), std::forward<Eps>(eps)
			);
		}

		template<class Scalar, class Function>
		auto make_diff_derivative(Function&& function, Scalar&& eps)
		{
			using Vec = Num::Arg::VecN<Scalar, 1>;
			using Mat = Num::Arg::MatNxM<Scalar, 1, 1>;

			auto wrapper = [&] (const Scalar& arg, const Vec& val)
			{
				return Vec(function(arg, val[0]));
			};

			return make_diff_jacobian<Vec, Mat>(wrapper, std::forward<Scalar>(eps));
		}
	}
}