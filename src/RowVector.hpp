#include "RowVector.h"

template<typename T>
RowVector<T>::RowVector() {

}

template<typename T>
RowVector<T>::~RowVector() {

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
inline void RowVector<T>::insert(size_t s, const T val)
{
	this->resize(s);
	for (size_t i = 0;i < s;++i) {
		getMatrix(*this).emplace_back(val);
	}
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
inline RowVector<T>& RowVector<T>::operator=(RowVector<T>&& m)
{
	Matrix<T>::operator=(std::move(m)); //no clue

	return *this;
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
inline RowVector<T>::RowVector(RowVector<T>&& c) : Matrix<T>(std::move(c))
{

}
