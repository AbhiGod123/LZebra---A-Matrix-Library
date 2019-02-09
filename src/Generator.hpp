#ifndef GENERATOR_HPP
#define GENERATOR_HPP
#include "Generator.h"

namespace gen {
	templ Matrix<T> eye(size_t row, size_t col) {
		Matrix<T> identity(row, col);

		for (size_t i = 0;i < row;++i) {
			identity(i, i) = 1;
		}

		return identity;
	}

	templ ColVector<T> unit(size_t size, size_t index)
	{
		ColVector<T> vec(size);
		vec(index) = 1;

		return vec;
	}
}

#endif // !GENERATOR_HPP