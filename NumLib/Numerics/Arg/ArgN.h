#pragma once


namespace Num
{
    namespace Arg
    {
        //// decls ////
        template<class Scalar, int N>
        class VecN;

        template<class Scalar, int N, int M>
        class MatNxM;

        template<class Scalar, int N>
        using MatNxN = MatNxM<Scalar, N, N>;



        //// VecN ////
        template<class ScalarType, int N>
        class VecN
        {
        public:
            using Scalar = ScalarType;


		public:
			static const int SIZE = N;


		public:
			template<class ... Init>
			VecN(const Init& ... initValues) : m_vec{initValues...}
			{
				static_assert(sizeof...(Init) == N, "Initialization values count should match vector size");
			}

            VecN() : m_vec()
            {}

            VecN(Scalar scalar)
            {
                for (int i = 0; i < N; i++)
                {
                    m_vec[i] = scalar;
                }
            }

            VecN(Scalar vec[N])
            {
                for (int i = 0; i < N; i++)
                {
                    m_vec[i] = vec[i];
                }
            }

            VecN(const VecN& vec)
            {
                for (int i = 0; i < N; i++)
                {
                    m_vec[i] = vec[i];
                }
            }

            VecN(std::initializer_list<Scalar> init)
            {
                int  i = 0;
                auto j = init.begin();

                while(i < N && j != init.end())
                {
                    m_vec[i] = *j;

                    i++;
                    j++;
                }

                while (i < N)
                {
                    m_vec[i] = Scalar(0);
                }
            }


            VecN operator + (const VecN& vec) const
            {
                VecN newVec;

                for (int i = 0; i < N; i++)
                {
                    newVec[i] = m_vec[i] + vec[i];
                }

                return newVec;
            }

            VecN operator - (const VecN& vec) const
            {
                VecN newVec;

                for (int i = 0; i < N; i++)
                {
                    newVec[i] = m_vec[i] - vec[i];
                }

                return newVec;
            }

            VecN operator * (const VecN& vec) const
            {
                VecN newVec;

                for (int i = 0; i < N; i++)
                {
                    newVec[i] = m_vec[i] * vec[i];
                }

                return newVec;
            }

            VecN operator / (const VecN& vec) const
            {
                VecN newVec;

                for (int i = 0; i < N; i++)
                {
                    newVec[i] = m_vec[i] / vec[i];
                }

                return newVec;
            }


            VecN operator + () const
            {
                return *this;
            }

            VecN operator - () const
            {
                VecN newVec;

                for (int i = 0; i < N; i++)
                {
                    newVec[i] = -m_vec[i];
                }

                return newVec;
            }


            VecN& operator = (const VecN& vec)
            {
                for (int i = 0; i < N; i++)
                {
                    m_vec[i] = vec[i];
                }

                return *this;
            }


            VecN& operator += (const VecN& vec)
            {
                for (int i = 0; i < N; i++)
                {
                    m_vec[i] += vec[i];
                }

                return *this;
            }

            VecN& operator -= (const VecN& vec)
            {
                for (int i = 0; i < N; i++)
                {
                    m_vec[i] -= vec[i];
                }

                return *this;
            }

            VecN& operator *= (const VecN& vec)
            {
                for (int i = 0; i < N; i++)
                {
                    m_vec[i] *= vec[i];
                }

                return *this;
            }

            VecN& operator /= (const VecN& vec)
            {
                for (int i = 0; i < N; i++)
                {
                    m_vec[i] /= vec[i];
                }

                return *this;
            }


            __forceinline Scalar& operator [] (int i)
            {
                return m_vec[i];
            }

            __forceinline const Scalar& operator [] (int i) const
            {
                return m_vec[i];
            }


            __forceinline const int size() const
            {
                return N;
            }


        protected:
            Scalar m_vec[N];
        };



        //// MatNxM ////
        template<class ScalarType, int N, int M>
        class MatNxM
        {
        public:
            using Scalar = ScalarType;

            template<class Type, int SIZE>
            using VectorType = VecN<Type, SIZE>;

			using Vector = VecN<Scalar, M>;


		public:
			static const int ROWS = N;
			static const int COLS = M;


		public:
            MatNxM() : m_mat()
            {}

            MatNxM(Scalar scalar)
            {
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < M; j++)
                    {
                        m_mat[i][j] = scalar;
                    }
                }
            }

            MatNxM(Scalar mat[N][M])
            {
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < M; j++)
                    {
                        m_mat[i][j] = mat[i][j];
                    }
                }
            }

            MatNxM(const MatNxM& mat)
            {
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < M; j++)
                    {
                        m_mat[i][j] = mat[i][j];
                    }
                }
            }

            MatNxM(std::initializer_list<Scalar> init)
            {
                int  i = 0;
                auto data = init.begin();

                while (i < N && data != init.end())
                {
                    int j = 0;

                    while (j < M && data != init.end())
                    {
                        m_mat[i][j] = *data;

                        j++;
                        data++;
                    }
                    while (j < M)
                    {
                        m_mat[i][j] = Scalar(0);
                    }

                    i++;
                }

                while (i < N)
                {
                    m_mat[i] = VecN<Scalar, M>();

                    i++;
                }
            }


            MatNxM operator + (const MatNxM& mat) const
            {
                MatNxM newMat;

                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < M; j++)
                    {
                        newMat[i][j] = m_mat[i][j] + mat[i][j];
                    }
                }

                return newMat;
            }

            MatNxM operator - (const MatNxM& mat) const
            {
                MatNxM newMat;

                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < M; j++)
                    {
                        newMat[i][j] = m_mat[i][j] - mat[i][j];
                    }
                }

                return newMat;
            }


            MatNxM operator + () const
            {
                return *this;
            }

            MatNxM operator - () const
            {
                MatNxM newMat;

                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < M; j++)
                    {
                        newMat[i][j] = -m_mat[i][j];
                    }
                }

                return newMat;
            }


            MatNxM& operator = (const MatNxM& mat)
            {
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < M; j++)
                    {
                        m_mat[i][j] = mat[i][j];
                    }
                }

                return *this;
            }


            MatNxM& operator += (const MatNxM& mat)
            {
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < M; j++)
                    {
                        m_mat[i][j] += mat[i][j];
                    }
                }

                return *this;
            }

            MatNxM& operator -= (const MatNxM& mat)
            {
                for (int i = 0; i < N; i++)
                {
                    for (int j = 0; j < M; j++)
                    {
                        m_mat[i][j] -= mat[i][j];
                    }
                }

                return *this;
            }


            __forceinline VecN<Scalar, M>& operator [] (int i)
            {
                return m_mat[i];
            }

            __forceinline const VecN<Scalar, M>& operator [] (int i) const
            {
                return m_mat[i];
            }


            __forceinline const int rows() const
            {
                return N;
            }

            __forceinline const int cols() const
            {
                return M;
            }


        protected:
            VecN<Scalar, M> m_mat[N];
        };



        //// operators : VecN & scalar ////
        template<class Scalar, int N>
        VecN<Scalar, N> operator + (const VecN<Scalar, N>& vec, Scalar scalar)
        {
            VecN<Scalar, N> newVec;

            for (int i = 0; i < N; i++)
            {
                newVec[i] = vec[i] + scalar;
            }

            return newVec;
        }

        template<class Scalar, int N>
        VecN<Scalar, N> operator - (const VecN<Scalar, N>& vec, Scalar scalar)
        {
            VecN<Scalar, N> newVec;

            for (int i = 0; i < N; i++)
            {
                newVec[i] = vec[i] - scalar;
            }

            return newVec;
        }

        template<class Scalar, int N>
        VecN<Scalar, N> operator * (const VecN<Scalar, N>& vec, Scalar scalar)
        {
            VecN<Scalar, N> newVec;

            for (int i = 0; i < N; i++)
            {
                newVec[i] = vec[i] * scalar;
            }

            return newVec;
        }

        template<class Scalar, int N>
        VecN<Scalar, N> operator / (const VecN<Scalar, N>& vec, Scalar scalar)
        {
            VecN<Scalar, N> newVec;

            for (int i = 0; i < N; i++)
            {
                newVec[i] = vec[i] / scalar;
            }

            return newVec;
        }


        template<class Scalar, int N>
        VecN<Scalar, N>& operator += (VecN<Scalar, N>& vec, Scalar scalar)
        {
            for (int i = 0; i < N; i++)
            {
                vec[i] += scalar;
            }

            return vec;
        }

        template<class Scalar, int N>
        VecN<Scalar, N>& operator -= (VecN<Scalar, N>& vec, Scalar scalar)
        {
            for (int i = 0; i < N; i++)
            {
                vec[i] -= scalar;
            }

            return vec;
        }

        template<class Scalar, int N>
        VecN<Scalar, N>& operator *= (VecN<Scalar, N>& vec, Scalar scalar)
        {
            for (int i = 0; i < N; i++)
            {
                vec[i] *= scalar;
            }

            return vec;
        }

        template<class Scalar, int N>
        VecN<Scalar, N>& operator /= (VecN<Scalar, N>& vec, Scalar scalar)
        {
            for (int i = 0; i < N; i++)
            {
                vec[i] /= scalar;
            }

            return vec;
        }


        //// operators : scalar & VecN ////
        template<class Scalar, int N>
        VecN<Scalar, N> operator * (Scalar scalar, const VecN<Scalar, N>& vec)
        {
            /*VecN<Scalar, N> newVec;

            for (int i = 0; i < N; i++)
            {
            newVec[i] = scalar * vec[i];
            }

            return newVec;*/
            return vec * scalar;
        }


        //// operators : MatNxM & scalars ////
        template<class Scalar, int N, int M>
        MatNxM<Scalar, N, M> operator + (const MatNxM<Scalar, N, N>& mat, Scalar scalar)
        {
            MatNxM<Scalar, N, M> newMat;

            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < M; j++)
                {
                    newMat[i][j] = mat[i][j] + scalar;
                }
            }

            return newMat;
        }

        template<class Scalar, int N, int M>
        MatNxM<Scalar, N, M> operator - (const MatNxM<Scalar, N, N>& mat, Scalar scalar)
        {
            MatNxM<Scalar, N, M> newMat;

            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < M; j++)
                {
                    newMat[i][j] = mat[i][j] - scalar;
                }
            }

            return newMat;
        }

        template<class Scalar, int N, int M>
        MatNxM<Scalar, N, M> operator * (const MatNxM<Scalar, N, N>& mat, Scalar scalar)
        {
            MatNxM<Scalar, N, M> newMat;

            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < M; j++)
                {
                    newMat[i][j] = mat[i][j] * scalar;
                }
            }

            return newMat;
        }

        template<class Scalar, int N, int M>
        MatNxM<Scalar, N, M> operator / (const MatNxM<Scalar, N, N>& mat, Scalar scalar)
        {
            MatNxM<Scalar, N, M> newMat;

            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < M; j++)
                {
                    newMat[i][j] = mat[i][j] / scalar;
                }
            }

            return newMat;
        }


        template<class Scalar, int N, int M>
        MatNxM<Scalar, N, M>& operator += (MatNxM<Scalar, N, N>& mat, Scalar scalar)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < M; j++)
                {
                    mat[i][j] += scalar;
                }
            }

            return mat;
        }

        template<class Scalar, int N, int M>
        MatNxM<Scalar, N, M>& operator -= (MatNxM<Scalar, N, N>& mat, Scalar scalar)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < M; j++)
                {
                    mat[i][j] -= scalar;
                }
            }

            return mat;
        }

        template<class Scalar, int N, int M>
        MatNxM<Scalar, N, M>& operator *= (MatNxM<Scalar, N, N>& mat, Scalar scalar)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < M; j++)
                {
                    mat[i][j] *= scalar;
                }
            }

            return mat;
        }

        template<class Scalar, int N, int M>
        MatNxM<Scalar, N, M>& operator /= (MatNxM<Scalar, N, N>& mat, Scalar scalar)
        {
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < M; j++)
                {
                    mat[i][j] /= scalar;
                }
            }

            return mat;
        }


        //// operators : scalar & MatNxM ////
        template<class Scalar, int N, int M>
        MatNxM<Scalar, N, M> operator * (Scalar scalar, const MatNxM<Scalar, N, M>& mat)
        {
            return mat * scalar;
        }


        //// operators : VecN & MatNxM ////
        template<class Scalar, int N, int M>
        VecN<Scalar, M> operator * (const VecN<Scalar, N>& vec, const MatNxM<Scalar, N, M>& mat)
        {
            VecN<Scalar, M> newVec;

            for (int j = 0; j < N; j++)
            {
                newVec += vec[j] * mat[j];
            }

            return newVec;
        }

        template<class Scalar, int N, int M>
        VecN<Scalar, N> operator * (const MatNxM<Scalar, N, M>& mat, const VecN<Scalar, M>& vec)
        {
            VecN<Scalar, N> newVec;

            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < M; j++)
                {
                    newVec[i] += mat[i][j] * vec[j];
                }
            }

            return newVec;
        }


        //// operators: MatNxM & MatNxM ////
        template<class Scalar, int N, int K, int M>
        MatNxM<Scalar, N, M> operator * (const MatNxM<Scalar, N, K>& mat1, const MatNxM<Scalar, K, M>& mat2)
        {
            MatNxM<Scalar, N, M> mat;

            for (int i = 0; i < N; i++)
            {
                mat[i] = mat1[i] * mat2;
            }

            return mat;
        }
    }
}