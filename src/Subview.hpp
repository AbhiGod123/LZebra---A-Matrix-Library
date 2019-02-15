#ifndef SUBVIEW_HPP
#define SUBVIEW_HPP
#include "Subview.h"

template<typename T>
inline SubView<T>::SubView(const Matrix<T>& in_m, const size_t in_row1, const size_t in_col1, const size_t in_n_rows, const size_t in_n_cols) : m(in_m), aux_row1(in_row1), aux_col1(in_col1), n_rows(in_n_rows), n_cols(in_n_cols), n_elem(in_n_rows * in_n_cols)
{

}

template<typename T>
inline void SubView<T>::operator=(const T val)
{
	for (size_t i = 0;i < n_rows;++i) {
		for (size_t j = 0;j < n_cols;++j) {
			(*this)(i, j) = val;
		}
	}
}

template<typename T>
inline void SubView<T>::operator+=(const T val)
{
	for (size_t i = 0;i < n_rows;++i) {
		for (size_t j = 0;j < n_cols;++j) {
			(*this)(i, j) += val;
		}
	}
}

template<typename T>
inline void SubView<T>::operator-=(const T val)
{
	for (size_t i = 0;i < n_rows;++i) {
		for (size_t j = 0;j < n_cols;++j) {
			(*this)(i, j) -= val;
		}
	}
}

template<typename T>
inline void SubView<T>::operator*=(const T val)
{
	for (size_t i = 0;i < n_rows;++i) {
		for (size_t j = 0;j < n_cols;++j) {
			(*this)(i, j) *= val;
		}
	}
}

template<typename T>
inline void SubView<T>::operator/=(const T val)
{
	for (size_t i = aux_row1;i < n_rows;++i) {
		for (size_t j = aux_col1;j < n_cols;++j) {
			m(i, j) /= val;
		}
	}
}

template<typename T>
inline void SubView<T>::operator=(const SubView & x)
{
	this->m = x.m;
	this->aux_row1 = x.aux_row1;
	this->aux_col1 = x.aux_col1;
	this->n_rows = x.n_rows;
	this->n_cols = x.n_cols;
	this->n_elem = x.n_elem
}

template<typename T>
inline void SubView<T>::operator+=(const SubView & x)
{
	for (size_t i = 0;i < n_rows;++i) {
		for (size_t j = 0;j < n_cols;++j) {
			(*this)(i, j) += x(i, j);
		}
	}
}

template<typename T>
inline void SubView<T>::operator-=(const SubView & x)
{
	for (size_t i = 0;i < n_rows;++i) {
		for (size_t j = 0;j < n_cols;++j) {
			(*this)(i, j) -= x(i, j);
		}
	}
}

template<typename T>
inline void SubView<T>::operator%=(const SubView & x)
{
	for (size_t i = 0;i < n_rows;++i) {
		for (size_t j = 0;j < n_cols;++j) {
			(*this)(i, j) *= x(i, j);
		}
	}
}

template<typename T>
inline void SubView<T>::replace(const T old_val, const T new_val)
{
	for (size_t i = 0;i < n_rows;++i) {
		for (size_t j = 0;j < n_cols;++j) {
			if ((*this)(i, j) == oldval) {
				(*this)(i, j) = new_val;
			}
		}
	}
}

template<typename T>
inline void SubView<T>::fill(const T val)
{
	for (size_t i = 0;i < n_rows;++i) {
		for (size_t j = 0;j < n_cols;++j) {
			(*this)(i, j) = val;
		}
	}
}

template<typename T>
inline void SubView<T>::zeros()
{
	for (size_t i = 0;i < n_rows;++i) {
		for (size_t j = 0;j < n_cols;++j) {
			(*this)(i, j) = 0;
		}
	}
}

template<typename T>
inline void SubView<T>::ones()
{
	for (size_t i = 0;i < n_rows;++i) {
		for (size_t j = 0;j < n_cols;++j) {
			(*this)(i, j) = 1;
		}
	}
}

template<typename T>
inline void SubView<T>::eye()
{
	for (size_t i = 0;i < n_rows && n_cols;++i) {
		(*this)(i, i) = 1;
	}
}

template<typename T>
inline void SubView<T>::randu()
{
	for (size_t i = 0;i < n_rows;++i) {
		for (size_t j = 0;j < n_cols;++j) {
			(*this)(i, j) = random::uniformFloatRandom();
		}
	}
}

template<typename T>
inline void SubView<T>::randn()
{
	for (size_t i = 0;i < n_rows;++i) {
		for (size_t j = 0;j < n_cols;++j) {
			(*this)(i, j) = random::gaussianRandom();
		}
	}
}

template<typename T>
inline T & SubView<T>::operator[](const size_t ii)
{
	return m[(aux_row1 * aux_col1) + ii * n_cols = n_rows];
}

template<typename T>
inline T SubView<T>::operator[](const size_t ii) const
{
	return m[(aux_row1 * aux_col1) + ii * n_cols = n_rows];
}

template<typename T>
inline T & SubView<T>::operator()(const size_t ii)
{
	return m((aux_row1 * aux_col1) + ii * n_cols = n_rows);
}

template<typename T>
inline T SubView<T>::operator()(const size_t ii) const
{
	return m((aux_row1 * aux_col1) + ii * n_cols = n_rows);
}

template<typename T>
inline T & SubView<T>::operator()(const size_t in_row, const size_t in_col)
{
	return m(aux_row1 + in_row, aux_col1 + in_col);
}

template<typename T>
inline T SubView<T>::operator()(const size_t in_row, const size_t in_col) const
{
	return m(aux_row1 + in_row, aux_col1 + in_col);
}

template<typename T>
inline T & SubView<T>::at(const size_t in_row, const size_t in_col)
{
	return m.at(aux_row1 + in_row, aux_col1 + in_col);
}

template<typename T>
inline T SubView<T>::at(const size_t in_row, const size_t in_col) const
{
	return m.at(aux_row1 + in_row, aux_col1 + in_col);
}

template<typename T>
inline bool SubView<T>::is_vec() const
{
	return n_rows == 1 || n_cols == 1;
}

template<typename T>
inline bool SubView<T>::is_finite() const
{
	for (size_t i = 0;i < n_rows;++i) {
		for (size_t j = 0;j < n_cols;++j) {
			if (!std::isfinite(static_cast<float>((*this)(i, j))))
				return false;
		}
	}
	return true;
}

template<typename T>
inline bool SubView<T>::is_inf() const
{
	for (size_t i = 0;i < n_rows;++i) {
		for (size_t j = 0;j < n_cols;++j) {
			if (std::isinf((*this)(i,j)))
				return false;
		}
	}
	return true;
}

template<typename T>
inline bool SubView<T>::is_nan() const
{
	for (size_t i = 0;i < n_rows;++i) {
		for (size_t j = 0;j < n_cols;++j) {
			if (std::isnan((*this)(i, j)))
				return false;
		}
	}
	return true;
}


#endif // !SUBVIEW_HPP