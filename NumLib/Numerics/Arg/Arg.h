#pragma once

#include <memory>


// decls ////
template<class Scalar>
class Vec;

template<class Scalar>
class Mat;



// Vec ////
template<class ScalarType>
class Vec
{
public:
    using UniquePtr = std::unique_ptr<Scalar[]>;
    using Scalar    = ScalarType;


    Vec() : m_vec(), m_size()
    {}

    Vec(int size) : Vec(size, Scalar())
    {}

    Vec(int size, Scalar scalar) : m_size(size), m_vec(new Scalar[size])
    {
        for (int i = 0; i < m_size; i++)
        {
            m_vec[i] = scalar;
        }
    }

    Vec(const Vec& vec) : m_size(vec.size()), m_vec(new Scalar[vec.size()])
    {
        for (int i = 0; i < m_size; i++)
        {
            m_vec[i] = vec[i];
        }
    }

    Vec(Vec&& vec)
    {
        m_size = vec.m_size;
        m_vec  = vec.m_vec;
    }


    Vec operator + (const Vec& vec) const
    {
        Vec newVec(*this);

        return newVec += vec;
    }

    Vec operator - (const Vec& vec) const
    {
        Vec newVec(*this);

        return newVec -= vec;
    }

    Vec operator * (const Vec& vec) const
    {
        Vec newVec(*this);

        return newVec *= vec;
    }

    Vec operator / (const Vec& vec) const
    {
        Vec newVec(*this);

        return newVec /= vec;
    }


    Vec operator + () const
    {
        return *this;
    }

    Vec operator - () const
    {
        Vec newVec(m_size);

        for (int i = 0; i < m_size; i++)
        {
            newVec[i] = -m_vec[i];
        }

        return newVec;
    }


    Vec& operator = (const Vec& vec)
    {
        if (this == &vec)
            return *this;

        m_size = vec.m_size;
        m_vec  = UniquePtr(new Scalar[m_size]);

        for (int i = 0; i < m_size; i++)
        {
            m_vec[i] = vec[i];
        }

        return *this;
    }


    Vec& operator += (const Vec& vec)
    {
        for (int i = 0; i < m_size; i++)
        {
            m_vec[i] += vec[i];
        }

        return *this;
    }

    Vec& operator -= (const Vec& vec)
    {
        for (int i = 0; i < m_size; i++)
        {
            m_vec[i] -= vec[i];
        }

        return *this;
    }

    Vec& operator *= (const Vec& vec)
    {
        for (int i = 0; i < m_size; i++)
        {
            m_vec[i] *= vec[i];
        }

        return *this;
    }

    Vec& operator /= (const Vec& vec)
    {
        for (int i = 0; i < m_size; i++)
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
        return m_size;
    }

protected:
    UniquePtr m_vec;
    int m_size;
};


// Mat ////
template<class ScalarType>
class Mat
{
public:
    using Scalar = ScalarType;


    Mat() : m_mat()
    {}

    Mat(int rows, int cols) : Mat(rows, cols, Scalar())
    {}

    Mat(int rows, int cols, Scalar scalar)
    {
        m_mat = Vec<Vec<Scalar>>(rows);
        for (int i = 0; i < m_rows; i++)
        {
            m_mat[i] = Vec<Scalar>(cols);

            for (int j = 0; j < m_cols; j++)
            {
                m_mat[i][j] = scalar;
            }
        }
    }

    Mat(const Mat& mat)
    {
        m_rows = mat.m_rows;
        m_cols = mat.m_cols;

        m_mat = Vec<Vec<Scalar>>(m_rows);
        for (int i = 0; i < m_rows; i++)
        {
            m_mat[i] = Vec<Scalar>(m_cols);

            for (int j = 0; j < m_cols; j++)
            {
                m_mat[i][j] = mat[i][j];
            }
        }
    }

    Mat(Mat&& mat)
    {
        m_mat  = mat.m_mat;
        m_cols = mat.m_cols;
        m_rows = mat.m_rows;
    }


    Mat operator + (const Mat& mat) const
    {
        Mat newMat(*this);

        return newMat += mat;
    }

    Mat operator - (const Mat& mat) const
    {
        Mat newMat(*this);

        return newMat -= mat;
    }

    Mat operator * (const Mat& mat) const
    {
        Mat<Scalar> newMat(m_rows, 0);

        for (int i = 0; i < m_rows; i++)
        {
            newMat[i] = m_mat[i] * mat;
        }

        return newMat;
    }


    Mat operator + () const
    {
        return *this;
    }

    Mat operator - () const
    {
        Mat newMat(m_rows, m_cols);

        for (int i = 0; i < m_rows; i++)
        {
            for (int j = 0; j < m_cols; j++)
            {
                newMat[i][j] = -m_mat[i][j];
            }
        }

        return newMat;
    }


    Mat& operator = (const Mat& mat)
    {
        m_rows = mat.m_rows;
        m_cols = mat.m_cols;
        
        m_mat = Vec<Vec<Scalar>>(m_rows);
        for (int i = 0; i < m_rows; i++)
        {
            m_mat[i] = Vec<Scalar>(m_cols);
            for (int j = 0; j < m_cols; j++)
            {
                m_mat[i][j] = mat[i][j];
            }
        }

        return *this;
    }


    Mat& operator += (const Mat& mat)
    {
        for (int i = 0; i < m_rows; i++)
        {
            for (int j = 0; j < m_cols; j++)
            {
                m_mat[i][j] += mat[i][j];
            }
        }

        return *this;
    }

    Mat& operator -= (const Mat& mat)
    {
        for (int i = 0; i < m_rows; i++)
        {
            for (int j = 0; j < m_cols; j++)
            {
                m_mat[i][j] -= mat[i][j];
            }
        }

        return *this;
    }


    __forceinline Vec<Scalar>& operator [] (int i)
    {
        return m_mat[i];
    }

    __forceinline const Vec<Scalar>& operator [] (int i) const
    {
        return m_mat[i];
    }


    __forceinline const int rows() const
    {
        return m_rows;
    }

    __forceinline const int cols() const
    {
        return m_cols;
    }


protected:
    Vec<Vec<Scalar>> m_mat;

    int m_rows;
    int m_cols;
};



// operators : Vec & scalar ////
template<class Scalar>
Vec<Scalar> operator + (const Vec<Scalar>& vec, Scalar scalar)
{
    Vec<Scalar> newVec(vec);

    for (int i = 0; i < vec.size(); i++)
    {
        newVec[i] += scalar;
    }

    return newVec;
}

template<class Scalar>
Vec<Scalar> operator - (const Vec<Scalar>& vec, Scalar scalar)
{
    Vec<Scalar> newVec(vec);

    for (int i = 0; i < vec.size(); i++)
    {
        newVec[i] -= scalar;
    }

    return newVec;
}

template<class Scalar>
Vec<Scalar> operator * (const Vec<Scalar>& vec, Scalar scalar)
{
    Vec<Scalar> newVec(vec);

    for (int i = 0; i < vec.size(); i++)
    {
        newVec[i] *= scalar;
    }

    return newVec;
}

template<class Scalar>
Vec<Scalar> operator / (const Vec<Scalar>& vec, Scalar scalar)
{
    Vec<Scalar> newVec(vec);

    for (int i = 0; i < vec.size(); i++)
    {
        newVec[i] /= scalar;
    }

    return newVec;
}


template<class Scalar>
Vec<Scalar>& operator += (Vec<Scalar>& vec, Scalar scalar)
{
    for (int i = 0; i < vec.size(); i++)
    {
        vec[i] += scalar;
    }

    return vec;
}

template<class Scalar>
Vec<Scalar>& operator -= (Vec<Scalar>& vec, Scalar scalar)
{
    for (int i = 0; i < vec.size(); i++)
    {
        vec[i] -= scalar;
    }

    return vec;
}

template<class Scalar>
Vec<Scalar>& operator *= (Vec<Scalar>& vec, Scalar scalar)
{
    for (int i = 0; i < vec.size(); i++)
    {
        vec[i] *= scalar;
    }

    return vec;
}

template<class Scalar>
Vec<Scalar>& operator /= (Vec<Scalar>& vec, Scalar scalar)
{
    for (int i = 0; i < vec.size(); i++)
    {
        vec[i] /= scalar;
    }

    return vec;
}


// operators : scalar & Vec ////
template<class Scalar>
Vec<Scalar> operator * (Scalar scalar, const Vec<Scalar>& vec)
{
    return vec * scalar;
}


// operators : Mat & scalars ////
template<class Scalar>
Mat<Scalar> operator + (const Mat<Scalar>& mat, Scalar scalar)
{
    Mat<Scalar> newMat(mat);

    return newMat += scalar;
}

template<class Scalar>
Mat<Scalar> operator - (const Mat<Scalar>& mat, Scalar scalar)
{
    Mat<Scalar> newMat(mat);

    return newMat -= scalar;
}

template<class Scalar>
Mat<Scalar> operator * (const Mat<Scalar>& mat, Scalar scalar)
{
    Mat<Scalar> newMat(mat);

    return newMat *= scalar;
}

template<class Scalar>
Mat<Scalar> operator / (const Mat<Scalar>& mat, Scalar scalar)
{
    Mat<Scalar> newMat(mat);

    return newMat /= scalar;
}


template<class Scalar>
Mat<Scalar>& operator += (Mat<Scalar>& mat, Scalar scalar)
{
    for (int i = 0; i < mat.rows(); i++)
    {
        mat[i] += scalar;
    }

    return mat;
}

template<class Scalar>
Mat<Scalar>& operator -= (Mat<Scalar>& mat, Scalar scalar)
{
    for (int i = 0; i < mat.rows(); i++)
    {
        mat[i] -= scalar;
    }

    return mat;
}

template<class Scalar>
Mat<Scalar>& operator *= (Mat<Scalar>& mat, Scalar scalar)
{
    for (int i = 0; i < mat.rows(); i++)
    {
        mat[i] *= scalar;
    }

    return mat;
}

template<class Scalar>
Mat<Scalar>& operator /= (Mat<Scalar>& mat, Scalar scalar)
{
    for (int i = 0; i < mat.rows(); i++)
    {
        mat[i] /= scalar;
    }

    return mat;
}


//// operators : scalar & Mat ////
template<class Scalar>
Mat<Scalar> operator * (Scalar scalar, const Mat<Scalar>& mat)
{
    return mat * scalar;
}


//// operators : Vec & Mat ////
template<class Scalar>
Vec<Scalar> operator * (const Vec<Scalar>& vec, const Mat<Scalar>& mat)
{
    Vec<Scalar> newVec(mat.cols());

    for (int i = 0; i < mat.cols(); i++)
    {
        for (int j = 0; j < mat.rows(); j++)
        {
            newVec[i] += vec[j] * mat[j][i];
        }
    }

    return newVec;
}

template<class Scalar>
Vec<Scalar> operator * (const Mat<Scalar>& mat, const Vec<Scalar>& vec)
{
    Vec<Scalar> newVec(mat.rows());

    for (int i = 0; i < mat.rows(); i++)
    {
        for (int j = 0; j < mat.cols(); j++)
        {
            newVec[i] = mat[i][j] * vec[j];
        }
    }

    return newVec;
}
