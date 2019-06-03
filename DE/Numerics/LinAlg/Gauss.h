#pragma once


#include "Common.h"


namespace
{
    template<class Matrix>
    auto findAbsMax(Matrix& matrix, int start)
    {
        std::pair<int, int> res{start, start};

        auto max = matrix[start][start];

        for (int i = start; i < matrix.rows(); i++)
        {
            for (int j = start; j < matrix.cols(); j++)
            {
                if (matrix[i][j] > max)
                {
                    res = {i, j};
                    max = matrix[i][j];
                }
            }
        }

        return res;
    }


    template<class Term>
    void swapComponents(Term& term, int first, int second)
    {
        std::swap(term[first], term[second]);
    }


    template<class Matrix>
    void swapRows(Matrix& matrix, int first, int second)
    {
        std::swap(matrix[first], matrix[second]);
    }


    template<class Matrix>
    void swapCols(Matrix& matrix, int first, int second)
    {
        for (int i = 0; i < matrix.rows(); i++)
        {
            std::swap(matrix[i][first], matrix[i][second]);
        }
    }


    template<class Matrix, class Term>
    void step(Matrix& mat, Term& term, int stepNumber)
    {
        auto mainElem = mat[stepNumber][stepNumber];

        for (int i = stepNumber + 1; i < mat.rows(); i++)
        {
            auto coef = mat[i][stepNumber] / mainElem;

            for (int j = stepNumber; j < mat.cols(); j++)
            {
                mat[i][j] -= mat[stepNumber][j] * coef;
            }

            term[i] -= term[stepNumber] * coef;
        }
    }
}


namespace Num
{
    namespace LinAlg
    {
        template<
            int N
            , template<class Scalar, int SIZE> class VectorType = Arg::VecN
        >
        class Permutation
        {
        public:
            using Vector = VectorType<int, N>;


        public:
            Permutation()
            {
                reset();
            }


            void permute(int first, int second)
            {
                std::swap(m_permutation[first], m_permutation[second]);

                m_parity ^= 1;
            }

            void reset()
            {
                for (int i = 0; i < m_permutation.size(); i++)
                {
                    m_permutation[i] = i;
                }

                m_parity = 0;
            }


            const int parity() const
            {
                return m_parity;
            }

            const int size() const
            {
                return m_permutation.size();
            }


            auto operator[] (int i)
            {
                return m_permutation[i];
            }

            const auto operator[] (int i) const
            {
                return m_permutation[i];
            }


        private:
            Vector m_permutation;
            int m_parity;
        };


        template<
            class ScalarType
            , int N
            , template<class Scalar, int ROWS, int COLS> class MatrixType = Arg::MatNxM
            , template<class Scalar, int SIZE> class VectorType = Arg::VecN
        >
        class GaussEliminationSingle
        {
        public:
            using Scalar = ScalarType;
            using Matrix = MatrixType<Scalar, N, N>;
            using Vector = VectorType<Scalar, N>;
            

        public:
            GaussEliminationSingle()
                : m_colsPerm()
                , m_rowsPerm()
                , m_determinant()
            {}


        private:
            bool forward(Matrix& mat, Vector& term)
            {
                for (int i = 0; i < mat.rows(); i++)
                {
                    auto index = findAbsMax(mat, i);

                    if (abs(mat[index.first][index.second]) < EPS<Scalar>)
                    {
                        return false;
                    }

                    if (index.first != i)
                    {
                        swapRows(mat, i, index.first);

                        swapComponents(term, i, index.first);

                        m_rowsPerm.permute(i, index.first);
                    }

                    if (index.second != i)
                    {
                        swapCols(mat, i, index.second);

                        m_colsPerm.permute(i, index.second);
                    }

                    step(mat, term, i);
                }

                return true;
            }

            void backward(Matrix& mat, Vector& term)
            {
                for (int i = mat.rows() - 1; i >= 0; i--)
                {
                    term[i] /= mat[i][i];

                    for (int j = i + 1; j < mat.cols(); j++)
                    {
                        term[i] -= term[j] * (mat[i][j] / mat[i][i]);
                    }
                }
            }

            void restoreComponentOrder(Vector& solution)
            {
                for (int i = 0; i < m_colsPerm.size(); )
                {
                    const int num = m_colsPerm[i];

                    if (num != i)
                    {
                        swapComponents(solution, i, num);
                        m_colsPerm.permute(i, num);
                    }
                    else
                    {
                        i++;
                    }
                }
            }

            void computeDeterminant(Matrix& mat)
            {
                for (int i = 0; i < mat.rows(); i++)
                {
                    m_determinant *= mat[i][i];
                }

                m_determinant *= ((m_colsPerm.parity() + m_rowsPerm.parity()) % 2 == 0 ? Scalar(1) : Scalar(-1));
            }


        public:
            void reset()
            {
                m_colsPerm.reset();
                m_rowsPerm.reset();

                m_determinant = Scalar(1);
            }


            bool solve(Matrix& mat, Vector& term)
            {
                if (!forward(mat, term))
                {
                    m_determinant = Scalar(0);

                    return false;
                }
                backward(mat, term);

                computeDeterminant(mat);

                restoreComponentOrder(term);

				reset();

                return true;
            }


            const Scalar determinant()
            {
                return m_determinant;
            }


        private:
            Permutation<N, VectorType> m_colsPerm;
            Permutation<N, VectorType> m_rowsPerm;

            Scalar m_determinant;
        };

        
        template<
            class ScalarType
            , int N, int M
            , template<class Scalar, int ROWS, int COLS> class MatrixType = Arg::MatNxM
            , template<class Scalar, int SIZE> class VectorType = Arg::VecN
        >
        class GaussEliminationMultiple
        {
        public:
            using Scalar = ScalarType;
            using Matrix = MatrixType<Scalar, N, N>;
            using Term   = MatrixType<Scalar, N, M>;


        public:
            GaussEliminationMultiple()
                : m_colsPerm()
                , m_rowsPerm()
                , m_determinant()
            {}


        private:
            bool forward(Matrix& mat, Term& term)
            {
                for (int i = 0; i < mat.rows(); i++)
                {
                    auto index = findAbsMax(mat, i);

                    if (abs(mat[index.first][index.second]) < EPS<Scalar>)
                    {
                        return false;
                    }

                    if (index.first != i)
                    {
                        swapRows(mat, i, index.first);

                        swapComponents(term, i, index.first);

                        m_rowsPerm.permute(i, index.first);
                    }

                    if (index.second != i)
                    {
                        swapCols(mat, i, index.second);

                        m_colsPerm.permute(i, index.second);
                    }

                    step(mat, term, i);
                }

                return true;
            }

            void backward(Matrix& mat, Term& term)
            {
                for (int i = mat.rows() - 1; i >= 0; i--)
                {
                    term[i] /= mat[i][i];

                    for (int j = i + 1; j < mat.cols(); j++)
                    {
                        term[i] -= term[j] * (mat[i][j] / mat[i][i]);
                    }
                }
            }

            void restoreComponentOrder(Term& solution)
            {
                for (int i = 0; i < m_colsPerm.size(); )
                {
                    const int num = m_colsPerm[i];

                    if (num != i)
                    {
                        swapComponents(solution, i, num);
                        m_colsPerm.permute(i, num);
                    }
                    else
                    {
                        i++;
                    }
                }
            }

            void computeDeterminant(Matrix& mat)
            {
                for (int i = 0; i < mat.rows(); i++)
                {
                    m_determinant *= mat[i][i];
                }

                m_determinant *= ((m_colsPerm.parity() + m_rowsPerm.parity()) % 2 == 0 ? Scalar(1) : Scalar(-1));
            }


        public:
            void reset()
            {
                m_colsPerm.reset();
                m_rowsPerm.reset();

                m_determinant = Scalar(1);
            }


            bool solve(Matrix& mat, Term& term)
            {
                if (!forward(mat, term))
                {
                    m_determinant = Scalar(0);

                    return false;
                }
                backward(mat, term);

                computeDeterminant(mat);

                restoreComponentOrder(term);

				reset();

                return true;
            }


            const Scalar determinant() const
            {
                return m_determinant;
            }


        private:
            Permutation<N, VectorType> m_colsPerm;
            Permutation<N, VectorType> m_rowsPerm;

            Scalar m_determinant;
        };
    }
}