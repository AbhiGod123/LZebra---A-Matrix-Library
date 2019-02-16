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
			(*this)(i, j) /= val;
		}
	}
}

template<typename T>
inline void SubView<T>::operator=(const SubView & x)
{
	for (size_t i = aux_row1;i < n_rows;++i) {
		for (size_t j = aux_col1;j < n_cols;++j) {
			(*this)(i, j) /= x.m(i,j);
		}
	}
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
			if ((*this)(i, j) == old_val) {
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
			std::cout << (*this)(i, j)<<std::endl;
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
	const size_t in_col = ii % n_cols;
	const size_t in_row = ii / n_cols;

	const size_t index = (in_row + aux_row1)*m.getCols() + aux_col1 + in_col;

	return const_cast<Matrix<T>&>(m)[index];
}

template<typename T>
inline T SubView<T>::operator[](const size_t ii) const
{
	const size_t in_col = ii % n_cols;
	const size_t in_row = ii / n_cols;

	const size_t index = (in_row + aux_row1)*m.getCols() + aux_col1 + in_col;

	return m[index];
}

template<typename T>
inline T & SubView<T>::operator()(const size_t ii)
{
	const size_t in_col = ii % n_cols;
	const size_t in_row = ii / n_cols;

	const size_t index = (in_row + aux_row1)*m.getCols() + aux_col1 + in_col;

	return const_cast<Matrix<T>&>(m)(index);
}

template<typename T>
inline T SubView<T>::operator()(const size_t ii) const
{
	const size_t in_col = ii % n_cols;
	const size_t in_row = ii / n_cols;

	const size_t index = (in_row + aux_row1)*m.getCols() + aux_col1 + in_col;

	return m(index);
}

template<typename T>
inline T & SubView<T>::operator()(const size_t in_row, const size_t in_col)
{
	const size_t index = (in_row + aux_row1)*m.getCols() + aux_col1 + in_col;

	return const_cast<Matrix<T>&>(m)(index);
}

template<typename T>
inline T SubView<T>::operator()(const size_t in_row, const size_t in_col) const
{
	const size_t index = (in_row + aux_row1)*m.getCols() + aux_col1 + in_col;

	return m(index);
}

template<typename T>
inline T & SubView<T>::at(const size_t in_row, const size_t in_col)
{
	const size_t index = (in_row + aux_row1)*m.getCols() + aux_col1 + in_col;

	if (index >= n_elem) {
		std::cout << "Index out of bounds!" << '\n';
	}

	return const_cast<Matrix<T>&>(m).at(index);
}

template<typename T>
inline T SubView<T>::at(const size_t in_row, const size_t in_col) const
{
	const size_t index = (in_row + aux_row1)*m.getCols() + aux_col1 + in_col;

	if (index >= n_elem) {
		std::cout << "Index out of bounds!" << '\n';
	}

	return m.at(index);
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

template<typename T>
inline SubViewRow<T> SubView<T>::row(const size_t row_num)
{
	return SubView<T>(m,row_num,n_cols,1,n_cols);
}

template<typename T>
inline const SubViewRow<T> SubView<T>::row(const size_t row_num) const
{
	return SubViewRow<T>(m, row_num, n_cols, 1, n_cols);
}

template<typename T>
inline SubViewCol<T> SubView<T>::col(const size_t col_num)
{
	return SubView<T>(m,n_rows,col_num,n_rows,1);
}

template<typename T>
inline const SubViewCol<T> SubView<T>::col(const size_t col_num) const
{
	return SubViewCol<T>(m, n_rows, col_num, n_rows, 1);
}

#endif // !SUBVIEW_HPP