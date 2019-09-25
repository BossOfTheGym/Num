#ifndef IVP_UTILS_H 
#define IVP_UTILS_H

#include <utility>
#include <type_traits>


#include "../Equ/Neuton.h"
#include "../Equ/Utils.h"

namespace Num
{
	namespace Ivp
	{

		template<class Argument, class Value, class Function>
		class DiffDerivative : public Function
		{
		public:
			template<class FuncType>
			DiffDerivative(FuncType&& function, const Value& eps) 
				: Function(std::forward<FuncType>(function))
				, mEps(eps)
			{}

			DiffDerivative(const DiffDerivative& deriv) = default;
			DiffDerivative(DiffDerivative&& deriv)      = default;


		public:
			Value operator()(const Argument& arg, const Value& val)
			{
				auto bound = [&] (const Value& value)
				{
					return Function::operator()(arg, value);
				};

				return Equ::make_diff_derivative(bound, mEps)(val);
			}

		private:
			Value mEps;
		};

		template<class Argument, class Value, class Function>
		auto make_diff_derivative(Function&& function, Value&& eps)
		{
			return DiffDerivative<Argument, Value, std::remove_reference_t<std::remove_cv_t<Function>>>(
				std::forward<Function>(function), std::forward<Value>(eps)
			);
		}



		template< class Vector, class Matrix, class Function>
		class DiffJacobian : public Function
		{
			static_assert(Vector::SIZE == Matrix::COLS && Vector::SIZE == Matrix::ROWS, "Non-matching sizes");


		public:
			using Scalar = typename Vector::Scalar;


		public:
			template<class FuncType, class EpsType>
			DiffJacobian(FuncType&& function, EpsType& eps) 
				: Function(std::forward<FuncType>(function))
				, mEps(std::forward<EpsType>(eps))
			{}

			DiffJacobian(const DiffJacobian& dj) = default;
			DiffJacobian(DiffJacobian&& dj)      = default;


			Matrix operator() (const Scalar& arg, const Vector& val)
			{
				const int N = Vector::SIZE;

				auto bound = [&] (const Vector& value) -> decltype(auto)
				{
					return Function::operator()(arg, value);
				};

				return Equ::make_diff_jacobian<Vector, Matrix>(bound, mEps)(val);
			}


		private:
			Scalar mEps;
		};

		template<class Vector, class Matrix, class Eps, class Function>
		auto make_diff_jacobian(Function&& function, Eps&& eps)
		{
			return DiffJacobian<Vector, Matrix, std::remove_reference_t<std::remove_cv_t<Function>>>(
				std::forward<Function>(function), std::forward<Eps>(eps)
			);
		}

	}
}

#endif