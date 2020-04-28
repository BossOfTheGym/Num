#pragma once

#include <vector>
#include <algorithm>

namespace Num::spline
{
    template<class Float>
    class CubicSpline1D
    {
    public:
        class SplineNode
        {
        public:
            SplineNode(Float a = {}, Float b = {}, Float c = {}, Float d = {}, Float left = {}, Float right = {})
                : m_a{a}, m_b{b}, m_c{c}, m_d{d}, m_left{left}, m_right{right}
            {}

            SplineNode(const SplineNode&) = default;
            SplineNode(SplineNode&&)      = default;

            ~SplineNode() = default;

            SplineNode& operator = (const SplineNode&) = default;
            SplineNode& operator = (SplineNode&&)      = default;


        public:
            Float computeValue(Float arg)
            {
                Float dlv = m_left - arg;
                Float drv = m_right - arg;

                return m_a * drv * drv * drv - m_b * dlv * dlv * dlv - m_c * dlv + m_d * drv;
            }

            Float computeDerivative(Float arg)
            {
                Float dlv = m_left - arg;
                Float drv = m_right - arg;

                return  -3 * m_a * drv * drv + 3 * m_b * dlv * dlv + m_c - m_d;
            }


        public:
            Float getA() const
            {
                return m_a;
            }

            Float getB() const
            {
                return m_b;
            }

            Float getC() const
            {
                return m_c;
            }

            Float getD() const
            {
                return m_d;
            }

            Float getLeft() const
            {
                return m_left;
            }

            Float getRight() const
            {
                return m_right;
            }


            void setA(Float a)
            {
                m_a = a;
            }

            void setB(Float b)
            {
                m_b = b;
            }

            void setC(Float c)
            {
                m_c = c;
            }

            void setD(Float d)
            {
                m_d = d;
            }

            void setLeft(Float left)
            {
                m_left = left;
            }

            void setRight(Float right)
            {
                m_right = right;
            }


        private:
            Float m_a{};
            Float m_b{};
            Float m_c{};
            Float m_d{};

            Float m_left{};
            Float m_right{};
        };


    private:
        struct SystemLine
        {
            Float a{};
            Float b{};
            Float c{};
            Float d{};
        };


    public:
        CubicSpline1D(const std::vector<Float>& args, const std::vector<Float>& values, Float firstMoment = static_cast<Float>(0), Float lastMoment = static_cast<Float>(0))
        {
            const int count   = std::min(args.size(), values.size());
            const int last    = (int)count - 1;
            const int prelast = (int)count - 2;

            // initialize systme
            std::vector<SystemLine> moments(count);

            // boundary condition
            moments[0].b = static_cast<Float>(1);    
            moments[0].d = firstMoment; 

            moments.back().b = static_cast<Float>(1);
            moments.back().d = lastMoment; 

            // initializing system
            for (int i = 1; i < last; i++)
            {
                auto hCurr = args[  i  ] - args[i - 1];
                auto hNext = args[i + 1] - args[  i  ];

                moments[i].a = hCurr / 6;
                moments[i].b = (hNext + hCurr) / 3;
                moments[i].c = hNext / 6;
                moments[i].d = (values[i + 1] - values[i]) / hNext - (values[i] - values[i - 1]) / hCurr;
            }

            // finding solution of a linear system
            // forward pass
            for (int i = 0; i < last; i++)
            {
                auto b = moments[i].b;
                moments[i].b = static_cast<Float>(1);
                moments[i].c /= b;
                moments[i].d /= b;

                auto a = moments[i + 1].a;
                moments[i + 1].a = static_cast<Float>(0);
                moments[i + 1].b -= moments[i].c * a;
                moments[i + 1].d -= moments[i].d * a;
            }
            { //last
                auto b = moments.back().b;
                moments.back().b = static_cast<Float>(1);
                moments.back().d /= b;
            }

            // backward pass
            for (int i = last; i > 1; i--)
            {
                auto c = moments[i - 1].c;
                moments[i - 1].c = static_cast<Float>(0);
                moments[i - 1].d -= moments[i].d * c;
            }

            //alloc
            m_splineNodes.resize(count - 1);

            //computing spline polynomial coefs
            for (int i = 0; i < count - 1; i++)
            {
                auto mCurr = moments[  i  ].d;
                auto mNext = moments[i + 1].d;
                auto h = args[i + 1] - args[  i  ];

                auto a = mCurr / 6 / h;
                auto b = mNext / 6 / h;
                auto c = (values[i + 1] - mNext * h * h / 6) / h;
                auto d = (values[  i  ] - mCurr * h * h / 6) / h;

                auto left  = args[  i  ];
                auto right = args[i + 1];

                m_splineNodes[i] = SplineNode(a, b, c, d, left, right);
            }

            m_left  = args.front();
            m_right = args.back();
            m_valueLeft  = values.front();
            m_valueRight = values.back();
        }


    public:
        Float operator() (Float arg)
        {
            if (arg <= m_left)
            {
                return m_valueLeft;
            }
            if (arg >= m_right)
            {
                return m_valueRight;
            }

            auto it = std::lower_bound(m_splineNodes.begin(), m_splineNodes.end(), arg,
                [] (auto& splineNode, auto& arg)
                {
                    return splineNode.getRight() < arg;
                }
            );

            return it->computeValue(arg);
        }


    private:
        std::vector<SplineNode> m_splineNodes{};

        Float m_left{};
        Float m_right{};
        Float m_valueLeft{};
        Float m_valueRight{};
    };
}