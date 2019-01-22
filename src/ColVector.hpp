#include "ColVector.h"

template<typename T>
ColVector<T>::ColVector() {

}

template<typename T>
ColVector<T>::~ColVector() {

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
inline void ColVector<T>::insert(size_t s, const T val)
{
	this->resize(s);
	for (size_t i = 0;i < s;++i) {
		getMatrix(*this).emplace_back(val);
	}
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
inline ColVector<T>& ColVector<T>::operator=(ColVector<T>&& m)
{
	Matrix<T>::operator=(std::move(m)); //no clue

	return *this;
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
inline ColVector<T>::ColVector(ColVector<T>&& c) : Matrix<T>(std::move(c))
{

}
