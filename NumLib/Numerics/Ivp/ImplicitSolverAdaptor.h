#pragma once

#include "Utils.h"

namespace Num
{
	namespace Ivp
	{
		template<class Argument, class Value, class ImplicitSolver>
		class ImplicitSolverAdaptor
		{
		public:
			using JacobianVecType = typename ImplicitSolver::JacobianVecType;
			using JacobianMatType = typename ImplicitSolver::JacobianMatType;

		public:
			template<class ImplicitSolverT>
			ImplicitSolverAdaptor(ImplicitSolverT solver)
				: m_solver(std::forward<ImplicitSolverT>(solver))
			{}

			ImplicitSolverAdaptor(const ImplicitSolver& solver)
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
		auto make_implicit_adaptor(ImplicitSolver&& solver)
		{
			using ImplicitSolverType = std::remove_reference_t<std::remove_cv_t<ImplicitSolver>>;

			return ImplicitSolverAdaptor<Argument, Value, ImplicitSolverType>(std::forward<ImplicitSolver>(solver));
		}
	}
}