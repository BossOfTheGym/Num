#ifndef Rect_h
#define Rect_h

#include <vector>
#include <cstdint>

namespace Utility
{
	using Size = size_t;

	template<class Scalar>
	class Rect : public std::vector<Scalar>
	{
	public:
		using Base = std::vector<Scalar>;


	private:
		template<class Element>
		class View
		{
		public: 
			using Pointer = Element*;
			using NoConst = std::remove_const_t<Element>;

		public:
			View(Pointer pointer, Size size) : m_pointer(pointer), m_size(size)
			{}

			Element& operator[] (Size i)
			{
				return m_pointer[i];
			}

			const NoConst& operator[] (Size i) const
			{
				return m_pointer[i];
			}

			Size size() const
			{
				return size;
			}


		private:
			Pointer m_pointer;
			Size m_size;
		};


	public:
		Rect(Size dim1, Size dim2, Scalar fill) : Base(dim1 * dim2, fill), m_dim1(dim1), m_dim2(dim2)
		{}


		Size dimension1() const
		{
			return m_dim1;
		}

		Size dimension2() const
		{
			return m_dim2;
		}


		auto operator[] (Size i)
		{
			return View(Base::data() + i * m_dim2, m_dim2);
		}

		auto operator[] (Size i) const
		{
			return View(Base::data() + i * m_dim2, m_dim2);
		}


	private:
		Size m_dim1;
		Size m_dim2;
	};
}

#endif