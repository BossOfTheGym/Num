#ifndef IMPLICIT_SOLVER_ADAPTER_H
#define IMPLICIT_SOLVER_ADAPTER_H

#include "Utils.h"

namespace Num
{
	namespace Ivp
	{
		template<class Argument, class Value, class ImplicitSolver>
		class ImplicitSolverAdapter
		{
		public:
			using JacobianVecType = typename ImplicitSolver::JacobianVecType;
			using JacobianMatType = typename ImplicitSolver::JacobianMatType;

		public:
			ImplicitSolverAdapter(ImplicitSolver&& solver)
				: m_solver(std::move(solver))
			{}

			ImplicitSolverAdapter(const ImplicitSolver& solver)
				: m_solver(solver)
			{}


		public:
			template<class Function>
			Value solve(
				Function&& func
				, const Argument& arg0
				, const Value& val0
				, const Argument& h
			)
			{
				return m_solver.solve(
					std::forward<Function>(func)
					, Num::Ivp::make_diff_jacobian<JacobianVecType, JacobianMatType>(func, h)
					, arg0
					, val0
					, h
				);
			}


		private:
			ImplicitSolver m_solver;
		};

		template<class Argument, class Value, class ImplicitSolver>
		auto make_implicit_adapter(ImplicitSolver&& solver)
		{
			return ImplicitSolverAdapter<Argument, Value, ImplicitSolver>(std::forward<ImplicitSolver>(solver));
		}
	}
}

#endif