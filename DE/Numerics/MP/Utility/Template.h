#ifndef Template_h
#define Template_h

#include <array>

namespace Utility
{
	template<class Scalar, size_t Dim1, size_t Dim2>
	using Template = std::array<std::array<Scalar, Dim2>, Dim1>;
}

#endif