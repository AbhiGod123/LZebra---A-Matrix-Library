#include "ColVector.h"

#ifndef COL_VECTOR_HPP
#define COL_VECTOR_HPP

template<typename T>
ColVector<T>::ColVector() : Matrix<T>(){

}

template<typename T>
inline void ColVector<T>::set_size(const size_t size)
{
	Matrix<T>::set_size(size, 1);
}

template<typename T>
inline void ColVector<T>::reshape(const size_t size)
{
	Matrix<T>::reshape(size, 1);
}

template<typename T>
inline void ColVector<T>::resize(const size_t size)
{
	Matrix<T>::resize(size, 1);
}

template<typename T>
inline void ColVector<T>::insert(const ColVector<T>& object)
{
	auto& mat = Matrix<T>::matrix;
	mat.insert(mat.end(), object.matrix.begin(), object.matrix.end());

	this->reshape(object.size + Matrix<T>::size);
}

template<typename T>
inline void ColVector<T>::insert(size_t s, const T val)
{
	auto& mat = Matrix<T>::matrix;
	mat.insert(mat.end(), s, val);

	this->reshape(s + Matrix<T>::size);
}

template<typename T>
inline void ColVector<T>::insert(const T val)
{
	auto& mat = Matrix<T>::matrix;
	mat.push_back(val);

	this->reshape(1 + Matrix<T>::size);
}

template<typename T>
inline void ColVector<T>::insert_zeros(size_t s)
{
	this->insert(s, 0);
}

template<typename T>
inline void ColVector<T>::insert_ones(size_t s)
{
	this->insert(s, 1);
}

template<typename T>
inline SubViewCol<T> ColVector<T>::get_col(size_t r1, size_t r2)
{
	if (r1 > r2)
		std::swap(r1, r2);

	if (r1 >= Matrix<T>::size || r2 >= Matrix<T>::size) {
		std::cout << "Col out of range!" << '\n';
	}

	const size_t subview_n_rows = r2 - r1 + 1;

	return SubViewCol<T>(Matrix<T>::matrix, 0, r1, subview_n_rows);
}

template<typename T>
inline SubViewCol<T> ColVector<T>::tail_col(size_t r1)
{
	return this->get_col(Matrix<T>::size - r1, Matrix<T>::size - 1);
}

template<typename T>
inline SubViewCol<T> ColVector<T>::head_col(size_t r2)
{
	return this->get_col(0, r2-1);
}

template<typename T>
inline void ColVector<T>::operator=(const Matrix<T>& x)
{
	Matrix<T>::operator= (x.getCols() == 1 ? x : tenseopr::vectorise(x, 0));
}

template<typename T>
inline ColVector<T>::ColVector(size_t size) : Matrix<T>(size, 1) {

}

template<typename T>
inline ColVector<T>::ColVector(size_t size, T elem) : Matrix<T>(size, 1, elem) {

}

template<typename T>
inline ColVector<T>::ColVector(const std::initializer_list<T>& list) : Matrix<T>(list.size(), 1, list)
{

}

template<typename T>
inline ColVector<T>::ColVector(const std::vector<T>& list) : Matrix<T>(list.size(), 1, list)
{

}

template<typename T>
inline ColVector<T>::ColVector(size_t size, const T * list) : Matrix<T>(size, 1, list)
{

}

template<typename T>
inline ColVector<T>::ColVector(const Matrix<T>& mat) : Matrix<T>(mat.getCols() == 1 ? mat : tenseopr::vectorise(mat, 0))
{
	
}

#endif // !COL_VECTOR_HPP