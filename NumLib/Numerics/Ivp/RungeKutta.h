#pragma once


#include "Common.h"

#include "../Equ/Neuton.h"


namespace Num
{
    namespace Ivp
    {
        template<
            class ScalarType
            , int N
            , template<class Scalar, int ROWS, int COLS> class MatrixType = Arg::MatNxM
            , template<class Scalar, int SIZE> class VectorType = Arg::VecN
        >
        struct ButcherTableau
        {
			template<class scalar, int rows, int cols>
			using MatTemplate = MatrixType<scalar, rows, cols>;

			template<class scalar, int size>
			using VecTemplate = VectorType<scalar, size>;


            using Scalar = ScalarType;
            using Vector = VectorType<Scalar, N>;
            using Matrix = MatrixType<Scalar, N, N>;

            static const int ORDER = N;


            Matrix mat; // coef matrix
            Vector c;   // time coefs
            Vector b;   // final coefs
        };


		template<class Scalar>
        struct Methods
        {
            //+++++ explicit methods +++++

            //tableau 3x3
			//------------------
			//  0  |  0   0   0  
			// 1/3 | 1/3  0   0
			// 2/3 |  0  2/3  0
			//-----|------------
			//     | 1/4  0  3/4
			//s = 3
			static ButcherTableau<Scalar, 3> heun3();


            //tableau 4x4
			//----------------------
			//  0  |  0   0   0   0  
			// 1/2 | 1/2  0   0   0
			// 1/2 |  0  1/2  0   0
			//  1  |  0   0   1   0
			//-----|----------------
			//     | 1/6 1/3 1/3 1/6
			//s = 4
			static ButcherTableau<Scalar, 4> classic4();


			//tableau 4x4
			//-------------------------------------------------------------
			//  0        |  0            0            0           0  
			// 0.4       | 0.4           0            0           0
			// 0.45573725| 0.29697761   0.15875964    0           0
			//  1        | 0.21810040  -3.05096516   3.83286476   0
			//-----------|-------------------------------------------------
			//           | 0.17476028  -0.55148066   1.20553560  0.17118478
			//s = 4
			static ButcherTableau<Scalar, 4> ralston4();


			//tableau 4x4
			//----------------------
			//  0  |  0    0   0   0  
			// 1/3 | 1/3   0   0   0
			// 2/3 |-1/3   1   0   0
			//  1  |  1   -1   1   0
			//-----|----------------
			//     | 1/8 3/8 3/8 1/8
			//s = 4
			static ButcherTableau<Scalar, 4> three_eighths_rule4();


            //+++++ implicit methods +++++

			//tableau 1x1
			//----------
			//  1  |  1
			//-----|----
			//     |  1
			//s = 1
			static ButcherTableau<Scalar, 1> backwordEuler1();

            //tableau 1x1
			//----------
			// 1/2 | 1/2
			//-----|----
			//     |  1
			//s = 2
			static ButcherTableau<Scalar, 1> midpoint2();


            //tableau 2x2
			// c1 = 1/2 - 1/6 * sqrt(3)
			// c2 = 1/2 + 1/6 * sqrt(3)
			//
			// a12 = 1/4 - 1/6 * sqrt(3)
			// a21 = 1/4 + 1/6 * sqrt(3)
			//
			//--------------
			// c1 | 1/4 a12
			// c2 | a21 1/4
			//----|---------
			//    | 1/2 1/2
			//s = 4
			static ButcherTableau<Scalar, 2> gaussLegendre4();


            //tableau 3x3
			// c1 = 1/2 - 1/10 * sqrt(15)
			// c3 = 1/2 + 1/10 * sqrt(15)
			//
			// a12 = 2/9 - 1/15 * sqrt(15)
			// a21 = 2/9 + 1/15 * sqrt(15)
			//
			// a13 = 5/36 - 1/30 * sqrt(15)
			// a31 = 5/36 + 1/30 * sqrt(15)
			//
			// a23 = 5/36 - 1/24 * sqrt(15)
			// a32 = 5/36 + 1/24 * sqrt(15)
			//
			//--------------------
			// c1  | 5/36 a12 a13  
			// 1/2 | a21  2/9 a23
			// c3  | a31  a32 5/36
			//-----|---------------
			//     | 5/18 4/9 5/18
			static ButcherTableau<Scalar, 3> gaussLegendre6();
        };

		template<class Scalar>
		ButcherTableau<Scalar, 3> Methods<Scalar>::heun3()
		{
			ButcherTableau<Scalar, 3> tableau;

			tableau.mat[1][0] = static_cast<Scalar>(1.0 / 3.0);
			tableau.mat[2][1] = static_cast<Scalar>(2.0 / 3.0);

			tableau.c[1] = static_cast<Scalar>(1.0 / 3.0);
			tableau.c[2] = static_cast<Scalar>(2.0 / 3.0);

			tableau.b[0] = static_cast<Scalar>(1.0 / 4.0);
			tableau.b[2] = static_cast<Scalar>(3.0 / 4.0);

			return tableau;
		}
		
		template<class Scalar>
		ButcherTableau<Scalar, 4> Methods<Scalar>::classic4()
		{
			ButcherTableau<Scalar, 4> tableau;

			tableau.mat[1][0] = static_cast<Scalar>(1.0 / 2.0);
			tableau.mat[2][1] = static_cast<Scalar>(1.0 / 2.0);
			tableau.mat[3][2] = static_cast<Scalar>(1.0);

			tableau.c[1] = static_cast<Scalar>(1.0 / 2.0);
			tableau.c[2] = static_cast<Scalar>(1.0 / 2.0);
			tableau.c[3] = static_cast<Scalar>(1.0);

			tableau.b[0] = static_cast<Scalar>(1.0 / 6.0);
			tableau.b[1] = static_cast<Scalar>(1.0 / 3.0);
			tableau.b[2] = static_cast<Scalar>(1.0 / 3.0);
			tableau.b[3] = static_cast<Scalar>(1.0 / 6.0);

			return tableau;
		}

		template<class Scalar>
		ButcherTableau<Scalar, 4> Methods<Scalar>::ralston4()
		{
			ButcherTableau<Scalar, 4> tableau;

			tableau.mat[1][0] = static_cast<Scalar>(0.4);

			tableau.mat[2][0] = static_cast<Scalar>(0.29697761);
			tableau.mat[2][1] = static_cast<Scalar>(0.15875964);

			tableau.mat[3][0] = static_cast<Scalar>(0.21810040);
			tableau.mat[3][1] = static_cast<Scalar>(-3.05096516);
			tableau.mat[3][2] = static_cast<Scalar>(3.83286476);


			tableau.c[1] = static_cast<Scalar>(0.4);
			tableau.c[2] = static_cast<Scalar>(0.45573725);
			tableau.c[3] = static_cast<Scalar>(1.0);


			tableau.b[0] = static_cast<Scalar>(0.17476028);
			tableau.b[1] = static_cast<Scalar>(-0.55148066);
			tableau.b[2] = static_cast<Scalar>(1.20553560);
			tableau.b[3] = static_cast<Scalar>(0.17118478);

			return tableau;
		}

		template<class Scalar>
		ButcherTableau<Scalar, 4> Methods<Scalar>::three_eighths_rule4()
		{
			ButcherTableau<Scalar, 4> tableau;

			tableau.mat[1][0] = static_cast<Scalar>(1.0 / 3.0);
			tableau.mat[2][0] = static_cast<Scalar>(-1.0 / 3.0);
			tableau.mat[2][1] = static_cast<Scalar>(1.0);
			tableau.mat[3][0] = static_cast<Scalar>(1.0);
			tableau.mat[3][1] = static_cast<Scalar>(-1.0);
			tableau.mat[3][2] = static_cast<Scalar>(1.0);

			tableau.c[1] = static_cast<Scalar>(1.0 / 3.0);
			tableau.c[2] = static_cast<Scalar>(2.0 / 3.0);
			tableau.c[3] = static_cast<Scalar>(1.0);

			tableau.b[0] = static_cast<Scalar>(1.0 / 8.0);
			tableau.b[1] = static_cast<Scalar>(3.0 / 8.0);
			tableau.b[2] = static_cast<Scalar>(3.0 / 8.0);
			tableau.b[3] = static_cast<Scalar>(1.0 / 8.0);

			return tableau;
		}


		template<class Scalar>
		ButcherTableau<Scalar, 1> Methods<Scalar>::backwordEuler1()
		{
			ButcherTableau<Scalar, 1> tableau;

			tableau.mat[0][0] = static_cast<Scalar>(1.0);

			tableau.c[0] = static_cast<Scalar>(1.0);

			tableau.b[0] = static_cast<Scalar>(1.0);

			return tableau;
		}

		template<class Scalar>
		ButcherTableau<Scalar, 1> Methods<Scalar>::midpoint2()
		{
			ButcherTableau<Scalar, 1> tableau;

			tableau.mat[0][0] = static_cast<Scalar>(1.0 / 2.0);

			tableau.c[0] = static_cast<Scalar>(1.0 / 2.0);

			tableau.b[0] = static_cast<Scalar>(1.0);

			return tableau;
		}

		template<class Scalar>
		ButcherTableau<Scalar, 2> Methods<Scalar>::gaussLegendre4()
		{
			ButcherTableau<Scalar, 2> tableau;

			tableau.mat[0][0] = static_cast<Scalar>(1.0 / 4.0);
			tableau.mat[0][1] = static_cast<Scalar>(1.0 / 4.0 - 1.0 / 6.0 * sqrt(3.0));
			tableau.mat[1][0] = static_cast<Scalar>(1.0 / 4.0 + 1.0 / 6.0 * sqrt(3));
			tableau.mat[1][1] = static_cast<Scalar>(1.0 / 4.0);

			tableau.c[0] = static_cast<Scalar>(1.0 / 2.0 - 1.0 / 6.0 * sqrt(3.0));
			tableau.c[1] = static_cast<Scalar>(1.0 / 2.0 + 1.0 / 6.0 * sqrt(3.0));

			tableau.b[0] = static_cast<Scalar>(1.0 / 2.0);
			tableau.b[1] = static_cast<Scalar>(1.0 / 2.0);

			return tableau;
		}

		template<class Scalar>
		ButcherTableau<Scalar, 3> Methods<Scalar>::gaussLegendre6()
		{
			ButcherTableau<Scalar, 3> tableau;

			tableau.mat[0][0] = static_cast<Scalar>(5.0 / 36.0);
			tableau.mat[2][2] = static_cast<Scalar>(5.0 / 36.0);

			tableau.mat[0][1] = static_cast<Scalar>(2.0 / 9.0  - 1.0 / 15.0 * sqrt(15.0));
			tableau.mat[2][1] = static_cast<Scalar>(2.0 / 9.0  + 1.0 / 15.0 * sqrt(15.0));

			tableau.mat[0][2] = static_cast<Scalar>(5.0 / 36.0 - 1.0 / 30.0 * sqrt(15.0));
			tableau.mat[2][0] = static_cast<Scalar>(5.0 / 36.0 + 1.0 / 30.0 * sqrt(15.0));

			tableau.mat[1][0] = static_cast<Scalar>(5.0 / 36.0 + 1.0 / 24.0 * sqrt(15.0));
			tableau.mat[1][2] = static_cast<Scalar>(5.0 / 36.0 - 1.0 / 24.0 * sqrt(15.0));

			tableau.mat[1][1] = static_cast<Scalar>(2.0 / 9.0);


			tableau.c[0] = static_cast<Scalar>(1.0 / 2.0 - 1.0 / 10.0 * sqrt(15));
			tableau.c[1] = static_cast<Scalar>(1.0 / 2.0);
			tableau.c[2] = static_cast<Scalar>(1.0 / 2.0 + 1.0 / 10.0 * sqrt(15));


			tableau.b[0] = static_cast<Scalar>(5.0 / 18.0);
			tableau.b[1] = static_cast<Scalar>(4.0 / 9.0);
			tableau.b[2] = static_cast<Scalar>(5.0 / 18.0);

			return tableau;
		}



        //unified for both Scalar & System
        template<
			  int SYSTEM_ORDER
            , class Argument
            , class Value
			, class Tableau = ButcherTableau<Argument, 4>
			, template<class Scalar, int SIZE> class VectorType = Arg::VecN
        >
        class RungeKuttaExplicit
        {
        public:
            using NextNode = Node<Argument, Value>;
			using Vector   = VectorType<Value, Tableau::ORDER>;


			RungeKuttaExplicit(Tableau&& tableau) : m_tableau(std::forward<Tableau>(tableau))
			{}

			template<class Function>
			NextNode solve(
				  Function&& func
				, const Argument& arg0
				, const Value& val0
				, const Argument& h
			)
			{
				const auto&[mat, cVec, bVec] = m_tableau;

				Vector kVec = Vector();
				for (int i = 0; i < Tableau::ORDER; i++)
				{
					for (int j = 0; j < i; j++)
					{
						kVec[i] += kVec[j] * mat[i][j];
					}

					kVec[i] = func(arg0 + h * cVec[i], val0 + h * kVec[i]);
				}

				Value sum = Value();
				for (int i = 0; i < Tableau::ORDER; i++)
				{
					sum += bVec[i] * kVec[i];
				}

				return {arg0 + h, val0 + h * sum};
			}
		

		private:
			Tableau m_tableau;
        };



        //unified for both Scalar & System
		template<
			  int SYSTEM_ORDER
			, class Argument
			, class Value   
			, class Tableau = ButcherTableau<Argument, 3>
			, class Solver  = Equ::NeutonSystem<Argument, Tableau::ORDER * SYSTEM_ORDER>
			, template<class Scalar, int SIZE>           class VectorType = Arg::VecN
			, template<class Scalar, int COLS, int ROWS> class MatrixType = Arg::MatNxM
		>
        class RungeKuttaImplicit
        {
        public:
			using NextNode = Node<Argument, Value>;


			RungeKuttaImplicit(
				  Tableau&& tableau
				, Argument&& eps
				, int iterationsLimit
			) 
				: m_tableau(std::forward<Tableau>(tableau))
				, m_solver(iterationsLimit, eps)
			{}


			template<class Function, class Jacobian>
			NextNode solve(
				  Function&& func
				, Jacobian&& jacobian
				, const Argument& arg0
				, const Value&    val0
				, const Argument& h
			)
			{
				//tableau parameters
				const auto&[mat, cVec, bVec] = m_tableau;


				//adjusted function
				using KVec = VectorType<Argument, Tableau::ORDER * SYSTEM_ORDER>;
				using KMat = VectorType<   Value, Tableau::ORDER>;

				auto specialFunc = [&] (const KVec& coefsVec) -> KVec
				{
					KVec  resultVec = KVec();
					KMat& resultMat = reinterpret_cast<KMat&>(resultVec);

					const KMat&  coefsMat = reinterpret_cast<const KMat&>(coefsVec);

					for (int i = 0; i < Tableau::ORDER; i++)
					{
						Value sums = Value();

						for (int j = 0; j < Tableau::ORDER; j++)
						{
							sums += mat[i][j] * coefsMat[j];
						}

						resultMat[i] = func(arg0 + h * cVec[i], val0 + h * sums);
						
					}

					return (coefsVec - resultVec);
				};


				//dirty "hack" to reinterpet result as it was matrix 
				using JacobianResult = MatrixType<Argument, SYSTEM_ORDER, SYSTEM_ORDER>;

				//adjusted jacobian
				using SpecialJacobian = MatrixType<Argument, Tableau::ORDER * SYSTEM_ORDER, Tableau::ORDER * SYSTEM_ORDER>;

				auto specialJacobian = [&] (const KVec& coefsVec) -> SpecialJacobian
				{
					// Special Jacobian dim = (Tableau::ORDER * BLOCK_DIM) * (Tableau::ORDER * BLOCK_DIM)
					//         block columns
					//        ________________
					//block  | _____   _____  |
					//lines  ||block| |block| |
					//       ||     | |     | |
					//       ||_____| |_____| |
					//       |                |
					//       | _____   _____  |+		

					//       ||block| |block| |
					//       ||     | |     | |
					//       ||_____| |_____| |
					//       |________________|
					//
					//      Block(in block)
					//      |   c o l u m n s|
					//      |l				 |
					//      |i				 |
					//      |n				 |
					//      |e				 |
					//      |s				 |
					//      |________________|
					// BLOCK_DIM = SYSTEM_ORDER * SYSTEM_ORDER

					const KMat& coefsMat = reinterpret_cast<const KMat&>(coefsVec);


					SpecialJacobian result;

					for (int blockLine = 0; blockLine < Tableau::ORDER; blockLine++)
					{
						Value sums = Value();

						for (int j = 0; j < Tableau::ORDER; j++)
						{
							sums += mat[blockLine][j] * coefsMat[j];
						}


						auto jacobianResult  = jacobian(arg0 + h * cVec[blockLine], val0 + h * sums);
						auto& systemJacobian = reinterpret_cast<JacobianResult&>(jacobianResult);

						for (int inBlockLine = 0; inBlockLine < SYSTEM_ORDER; inBlockLine++)
						{
							for(int blockColumn = 0; blockColumn < Tableau::ORDER; blockColumn++)
							{
								for (int inBlockColumn = 0; inBlockColumn < SYSTEM_ORDER; inBlockColumn++)
								{
									result[SYSTEM_ORDER * blockLine + inBlockLine][SYSTEM_ORDER * blockColumn + inBlockColumn] 
										= -systemJacobian[inBlockLine][inBlockColumn] * h * mat[blockLine][blockColumn];
								}
							}

							result[SYSTEM_ORDER * blockLine + inBlockLine][SYSTEM_ORDER * blockLine + inBlockLine] 
								+= static_cast<Argument>(1.0);
						}
					}

					return result;
				};


				//initial value
				KVec  initialVec = KVec();
				KMat& initialMat = reinterpret_cast<KMat&>(initialVec);

				Value val = func(arg0, val0);
				for (int i = 0; i < Tableau::ORDER; i++)
				{
					initialMat[i] = val;
				}

				auto solution     = m_solver.solve(std::move(specialFunc), std::move(specialJacobian), std::move(initialVec));
				KMat& solutionMat = reinterpret_cast<KMat&>(solution.first);


				//computing result
				Value sum = Value();
				for(int i = 0; i < Tableau::ORDER; i++)
				{
					sum += bVec[i] * solutionMat[i];
				}


				return {arg0 + h, val0 + h * sum};
			}


		public:
			Solver&  solver()
			{
				return m_solver;
			}

			Tableau& tableau()
			{
				return m_tableau;
			}


		private:
			Solver   m_solver;
			Tableau m_tableau;
        };



		template<class Argument, class Value, class Function>
		class DiffDerivative : public Function
		{
		public:
			DiffDerivative(Function&& function, Value&& eps) 
				: Function(std::forward<Function>(function))
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

				return Equ::DiffDerivative(bound, mEps)(val);
			}

		private:
			Value mEps;
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


			Matrix operator() (const Scalar& arg, const Vector& val)
			{
				const int N = Vector::SIZE;

				auto bound = [&] (const Vector& value) -> decltype(auto)
				{
					return Function::operator()(arg, value);
				};

				return Equ::DiffJacobian<Vector, Matrix, decltype(bound)>(std::move(bound), std::move(mEps))(val);
			}


		private:
			Scalar mEps;
		};
    }
}