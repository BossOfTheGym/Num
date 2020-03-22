#pragma once


#include <functional>

#include "../Common/Common.h"
#include "../Arg/ArgN.h"


namespace Num
{
    namespace Utils
    {
        enum NormType : int
        {
              Cube 
            , Octo
            , Sphere
        };

        template<class Vector, NormType type>
        class Norm;

        template<class Vector>
        class Norm<Vector, Cube>
        {
        public:
            using Scalar = typename Vector::Scalar;


            Scalar operator ()(const Vector& vec)
            {
                Scalar value = abs(vec[0]);

                for (int i = 1; i < vec.size(); i++)
                {
                    Scalar next  = abs(vec[i]);
                    value = (value > next ? value : next);
                }

                return value;
            }
        };

        template<class Vector>
        class Norm<Vector, Octo>
        {
        public:
            using Scalar = typename Vector::Scalar;


            Scalar operator()(const Vector& vec)
            {
                Scalar value = Scalar(0);

                for (int i = 0; i < vec.size(); i++)
                {
                    value += abs(vec[i]);
                }

                return value;
            }
        };

        template<class Vector>
        class Norm<Vector, Sphere>
        {
        public:
            using Scalar = typename Vector::Scalar;


            Scalar operator()(const Vector& vec)
            {
                Scalar value = Scalar(0);

                for (int i = 0; i < vec.size(); i++)
                {
                    value += vec[i] * vec[i];
                }

                return std::sqrt(value);
            }
        };
    }
}