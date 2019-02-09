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
	templ ColVector<T> linspace(size_t start, size_t end, size_t size)
	{
		ColVector<T> vec(size);
	
		return vec;
	}

	templ ColVector<T> ones(size_t size)
	{
		return ColVector<T>(size, 1);
	}

	templ Matrix<T> ones(size_t rows, size_t cols)
	{
		return Matrix<T>(rows, cols, 1);
	}
	templdef ColVector<T> randi(size_t size)
	{
		ColVector<T> vec(size);
		vec.randi();

		return vec;
	}
	templdef Matrix<T> randi(size_t rows, size_t cols)
	{
		Matrix<T> mat(rows, cols);
		mat.randi();

		return mat;
	}
	ftempldec ColVector<T> randu(size_t size)
	{
		ColVector<T> vec(size);
		vec.randu();

		return vec;
	}
	ftempldec Matrix<T> randu(size_t rows, size_t cols)
	{
		Matrix<T> mat(rows, cols);
		mat.randu();

		return mat;
	}
	ftempldec ColVector<T> randn(size_t size)
	{
		ColVector<T> vec(size);
		vec.randn();

		return vec;
	}
	ftempldec Matrix<T> randn(size_t rows, size_t cols)
	{
		Matrix<T> mat(rows, cols);
		mat.randn();

		return mat;
	}
	ftempldec ColVector<T> regspace(size_t start, size_t end)
	{
		ColVector<T> vec(end);

		for (size_t i = 0;i < vec.getSize();++i) {
			vec(i) = start + i * 1;
		}

		return vec;
	}
	ftempldec ColVector<T> regspace(size_t start, size_t delta, size_t end)
	{
		ColVector<T> vec(end);

		for (size_t i = 0;i < vec.getSize();++i) {
			vec(i) = start + i * delta;
		}

		return vec;
	}
	templ Matrix<T> toeplitz(ccvec vec)
	{
		Matrix<T> mat(vec.getSize(),vec.getSize());

		for (size_t i = 0;i < mat.getRows();++i) {
			for (size_t j = 0;j < mat.getCols();++j) {
				mat(i, j) = vec(::abs(static_cast<long long>(i) - j));

			}
		}

		return mat;
	}
	templ ColVector<T> zeros(size_t size)
	{
		return ColVector<T>(size, 0);
	}
	templ Matrix<T> zeros(size_t rows, size_t cols)
	{
		return Matrix<T>(rows, cols, 0);
	}
}

#endif // !GENERATOR_HPP