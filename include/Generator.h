#ifndef GENERATOR_H
#define GENERATOR_H
#include "Include.h"

namespace gen {
	templ Matrix<T> eye(size_t row, size_t col);
	templ ColVector<T> unit(size_t size, size_t index);
}

#endif // !GENERATOR_H