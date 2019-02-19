#ifndef GENERATOR_H
#define GENERATOR_H
#include "Include.h"

namespace LZebra {
	namespace gen {
		templ Matrix<T> eye(size_t row, size_t col);
		templ ColVector<T> unit(size_t size, size_t index);
		templ ColVector<T> linspace(size_t start, size_t end, size_t size = 100);
		templ ColVector<T> ones(size_t size);
		templ Matrix<T> ones(size_t rows, size_t cols);
		itempldec ColVector<T> randi(size_t size, T min = 0, T max = std::numeric_limits<T>::max);
		itempldec Matrix<T> randi(size_t rows, size_t cols, T min = 0, T max = std::numeric_limits<T>::max);
		ftempldec ColVector<T> randu(size_t size);
		ftempldec Matrix<T> randu(size_t rows, size_t cols);
		ftempldec ColVector<T> randn(size_t size);
		ftempldec Matrix<T> randn(size_t rows, size_t cols);
		templ ColVector<T> regspace(size_t start);
		templ ColVector<T> regspace(size_t start, size_t end);
		templ ColVector<T> regspace(size_t start, size_t delta, size_t end);
		templ Matrix<T> toeplitz(ccvec vec);
		templ ColVector<T> zeros(size_t size);
		templ Matrix<T> zeros(size_t rows, size_t cols);
	}
};
#endif // !GENERATOR_H