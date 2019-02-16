#include "RowVector.h"

#ifndef ROW_VECTOR_HPP
#define ROW_VECTOR_HPP

template<typename T>
RowVector<T>::RowVector() : Matrix<T>(){

}

template<typename T>
inline void RowVector<T>::set_size(const size_t size)
{
	Matrix<T>::set_size(1, size);
}

template<typename T>
inline void RowVector<T>::reshape(const size_t size)
{
	Matrix<T>::reshape(1, size);
}

template<typename T>
inline void RowVector<T>::resize(const size_t size)
{
	Matrix<T>::resize(1, size);
}

template<typename T>
inline void RowVector<T>::insert(const RowVector<T>& object)
{
	auto& mat = Matrix<T>::matrix;
	mat.insert(mat.end(), object.matrix.begin(), object.matrix.end());

	this->reshape(object.size + Matrix<T>::size);
}

template<typename T>
inline void RowVector<T>::insert(size_t s, const T val)
{
	auto& mat = Matrix<T>::matrix;
	mat.insert(mat.end(), s, val);

	this->reshape(s + Matrix<T>::size);
}

template<typename T>
inline void RowVector<T>::insert(const T val)
{
	auto& mat = Matrix<T>::matrix;
	mat.push_back(val);

	this->reshape(1 + Matrix<T>::size);
}

template<typename T>
inline void RowVector<T>::insert_zeros(size_t s)
{
	this->insert(0, s);
}

template<typename T>
inline void RowVector<T>::insert_ones(size_t s)
{
	this->insert(1, s);
}

template<typename T>
inline SubViewRow<T> RowVector<T>::get_row(size_t c1, size_t c2)
{
	if (c1 > c2)
		std::swap(c1, c2);

	if (c1 >= Matrix<T>::size || c2 >= Matrix<T>::size) {
		std::cout << "Col out of range!" << '\n';
	}

	const size_t subview_n_cols = c2 - c1 + 1;

	return SubViewRow<T>(m, c1, 0, subview_n_cols);
}

template<typename T>
inline SubViewRow<T> RowVector<T>::tail_row(size_t c1)
{
	return this->get_row(Matrix<T>::size - 1, Matrix<T>::size - c1);
}

template<typename T>
inline SubViewRow<T> RowVector<T>::head_row(size_t c2)
{
	return this->get_row(c2 - 1, 0);
}

template<typename T>
inline void RowVector<T>::operator=(const Matrix<T>& x)
{
	Matrix<T>::operator= (x.getRows() == 1 ? x : tenseopr::vectorise(x, 1));
}

template<typename T>
inline RowVector<T>::RowVector(size_t size) : Matrix<T>(1, size) {

}

template<typename T>
inline RowVector<T>::RowVector(size_t size, T elem) : Matrix<T>(1, size, elem) {

}

template<typename T>
inline RowVector<T>::RowVector(const std::initializer_list<T>& list) : Matrix<T>(1, list.size(), list)
{

}

template<typename T>
inline RowVector<T>::RowVector(const std::vector<T>& list) : Matrix<T>(1, list.size(), list)
{

}

template<typename T>
inline RowVector<T>::RowVector(size_t size, const T * list) : Matrix<T>(1, size, list)
{

}

template<typename T>
inline RowVector<T>::RowVector(const Matrix<T>& mat) : Matrix<T>(mat.getRows() == 1 ? mat : tenseopr::vectorise(mat,1))
{

}

#endif // !ROW_VECTOR_HPP