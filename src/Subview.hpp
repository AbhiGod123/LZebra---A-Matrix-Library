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
	for (size_t i = 0;i < n_rows;++i) {
		for (size_t j = 0;j < n_cols;++j) {
			(*this)(i, j) /= val;
		}
	}
}

template<typename T>
inline void SubView<T>::operator=(const SubView & x)
{
	for (size_t i = 0;i < n_rows;++i) {
		for (size_t j = 0;j < n_cols;++j) {
			(*this)(i, j) = x(i,j);
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
inline void SubView<T>::operator/=(const SubView<T>& x)
{
	for (size_t i = 0;i < n_rows;++i) {
		for (size_t j = 0;j < n_cols;++j) {
			(*this)(i, j) /= x(i, j);
		}
	}
}

template<typename T>
inline void SubView<T>::operator++()
{
	(*this) += 1;
}

template<typename T>
inline void SubView<T>::operator++(int)
{
	(*this) += 1;
}

template<typename T>
inline void SubView<T>::operator--()
{
	(*this) -= 1;
}

template<typename T>
inline void SubView<T>::operator--(int)
{
	(*this) -= 1;
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
inline size_t SubView<T>::getRows() const
{
	return n_rows;
}

template<typename T>
inline size_t SubView<T>::getCols() const
{
	return n_cols;
}

template<typename T>
inline size_t SubView<T>::getSize() const
{
	return n_elem;
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
			if (std::isinf(static_cast<float>((*this)(i, j))))
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
			if (std::isnan(static_cast<float>((*this)(i, j))))
				return false;
		}
	}
	return true;
}

template<typename T>
inline SubViewRow<T> SubView<T>::row(const size_t row_num)
{
	return SubViewRow<T>(m, aux_col1, row_num + aux_row1, n_cols);
}

template<typename T>
inline const SubViewRow<T> SubView<T>::row(const size_t row_num) const
{
	return SubViewRow<T>(m, aux_col1,row_num + aux_row1, n_cols);
}

template<typename T>
inline SubViewCol<T> SubView<T>::col(const size_t col_num)
{
	return SubViewCol<T>(m,col_num + aux_col1,aux_row1,n_rows);
}

template<typename T>
inline const SubViewCol<T> SubView<T>::col(const size_t col_num) const
{
	return SubViewCol<T>(m, col_num + aux_col1, aux_row1, n_rows);
}

template<typename T>
inline SubView<T> SubView<T>::rows(size_t in_row1, size_t in_row2)
{
	if (in_row2 > in_row2)
		std::swap(in_row1, in_row2);

	if (in_row1 > n_rows || in_row2 > n_rows)
	{
		std::cout << "Row out of range!" << '\n';
	}

	std::cout << "fdsfs" << std::endl;

	return SubView<T>(m, in_row1 + aux_row1, aux_col1, in_row2 - in_row1 + 1, n_cols);
}

template<typename T>
inline const SubView<T> SubView<T>::rows(size_t in_row1, size_t in_row2) const
{
	if (in_row2 > in_row2)
		std::swap(in_row1, in_row2);

	if (in_row1 > n_rows || in_row2 > n_rows)
	{
		std::cout << "Row out of range!" << '\n';
	}

	return SubView<T>(m, in_row1 + aux_row1, aux_col1, in_row2 - in_row1 + 1, n_cols);
}

template<typename T>
inline SubView<T> SubView<T>::cols(size_t in_col1, size_t in_col2)
{
	if (in_col2 > in_col2)
		std::swap(in_col1, in_col2);

	if (in_col1 > n_cols || in_col2 > n_cols)
	{
		std::cout << "Row out of range!" << '\n';
	}

	return SubView<T>(m, aux_row1, in_col1 + aux_col1, n_rows, in_col2 - in_col1 + 1);
}

template<typename T>
inline const SubView<T> SubView<T>::cols(size_t in_col1, size_t in_col2) const
{
	if (in_col2 > in_col2)
		std::swap(in_col1, in_col2);

	if (in_col1 > n_cols || in_col2 > n_cols)
	{
		std::cout << "Row out of range!" << '\n';
	}

	return SubView<T>(m, aux_row1, in_col1 + aux_col1, n_rows, in_col2 - in_col1 + 1);
}

template<typename T>
inline SubView<T> SubView<T>::submat(size_t r1, size_t c1, size_t r2, size_t c2)
{
	if (r1 > r2)
		std::swap(r1, r2);
	if (c1 > c2)
		std::swap(c1, c2);

	if (r1 >= n_rows || r2>= n_rows) {
		std::cout << "Row out of range!" << '\n';
	}
	if (c1 >= n_cols || c2>= n_cols) {
		std::cout << "Col out of range!" << '\n';
	}

	return SubView<T>(m, r1, c1, r2 - r1 + 1, c2 - c1 + 1);
}

template<typename T>
inline const SubView<T> SubView<T>::submat(size_t r1, size_t c1, size_t r2, size_t c2) const
{
	if (r1 > r2)
		std::swap(r1, r2);
	if (c1 > c2)
		std::swap(c1, c2);

	if (r1 >= n_rows || r2 >= n_rows) {
		std::cout << "Row out of range!" << '\n';
	}
	if (c1 >= n_cols || c2 >= n_cols) {
		std::cout << "Col out of range!" << '\n';
	}

	return SubView<T>(m, r1, c1, r2 - r1 + 1, c2 - c1 + 1);
}

template<typename T>
inline void SubView<T>::print() const
{
	for (size_t i = 0;i < n_rows;++i)
	{
		for (size_t f = 0;f < n_cols;++f)
		{
			std::cout << (*this)(i, f) << ' ';
		}
		std::cout << '\n';
	}
}

template<typename T>
inline SubViewCol<T>::SubViewCol(const Matrix<T>& in_m, const size_t in_col1, const size_t in_row1, const size_t in_n_rows) : SubView<T>(in_m, in_row1, in_col1, in_n_rows, 1)
{

}

template<typename T>
inline SubViewCol<T> SubViewCol<T>::subvec(const size_t in_row1, const size_t in_row2)
{
	const size_t subview_n_rows = in_row2 - in_row1 + 1;

	const size_t base_row1 = SubView<T>::aux_row1 + in_row1;

	return SubViewCol<T>(SubView<T>::m, SubView<T>::aux_col1, base_row1, subview_n_rows);
}

template<typename T>
inline const SubViewCol<T> SubViewCol<T>::subvec(const size_t in_row1, const size_t in_row2) const
{
	const size_t subview_n_rows = in_row2 - in_row1 + 1;

	const size_t base_row1 = SubView<T>::aux_row1 + in_row1;

	return SubViewCol<T>(SubView<T>::m, SubView<T>::aux_col1, base_row1, subview_n_rows);
}

template<typename T>
inline SubViewCol<T> SubViewCol<T>::head(const size_t N)
{
	return SubViewCol<T>(SubView<T>::m, SubView<T>::aux_col1, SubView<T>::aux_row1, N);
}

template<typename T>
inline const SubViewCol<T> SubViewCol<T>::head(const size_t N) const
{
	return SubViewCol<T>(SubView<T>::m, SubView<T>::aux_col1, SubView<T>::aux_row1, N);
}

template<typename T>
inline SubViewCol<T> SubViewCol<T>::tail(const size_t N)
{
	const size_t start_row = SubView<T>::aux_row1 + SubView<T>::n_rows - N;

	return SubViewCol<T>(SubView<T>::m, SubView<T>::aux_col1, start_row, N);
}

template<typename T>
inline const SubViewCol<T> SubViewCol<T>::tail(const size_t N) const
{
	const size_t start_row = SubView<T>::aux_row1 + SubView<T>::n_rows - N;

	return SubViewCol<T>(SubView<T>::m, SubView<T>::aux_col1, start_row, N);
}

template<typename T>
inline SubViewCol<T>::SubViewCol(const Matrix<T>& in_m, const size_t in_col) : SubView<T>(in_m, 0, in_col, in_m.getRows(), 1)
{

}

template<typename T>
inline SubViewRow<T> SubViewRow<T>::subvec(const size_t in_col1, const size_t in_col2)
{
	const size_t subview_n_cols = in_col2 - in_col1 + 1;

	const size_t base_col1 = SubView<T>::aux_col1 + in_col1;

	return SubViewRow<T>(SubView<T>::m, SubView<T>::aux_row1, base_col1, subview_n_cols);
}

template<typename T>
inline const SubViewRow<T> SubViewRow<T>::subvec(const size_t in_col1, const size_t in_col2) const
{
	const size_t subview_n_cols = in_col2 - in_col1 + 1;

	const size_t base_col1 = SubView<T>::aux_col1 + in_col1;

	return SubViewRow<T>(SubView<T>::m, SubView<T>::aux_row1, base_col1, subview_n_cols);
}

template<typename T>
inline SubViewRow<T> SubViewRow<T>::head(const size_t N)
{
	return SubViewRow<T>(SubView<T>::m, SubView<T>::aux_row1, SubView<T>::aux_col1, N);
}

template<typename T>
inline const SubViewRow<T> SubViewRow<T>::head(const size_t N) const
{
	return SubViewRow<T>(SubView<T>::m, SubView<T>::aux_row1, SubView<T>::aux_col1, N);
}

template<typename T>
inline SubViewRow<T> SubViewRow<T>::tail(const size_t N)
{
	const size_t start_col = SubView<T>::aux_col1 + SubView<T>::n_cols - N;

	return SubViewRow<T>(SubView<T>::m, SubView<T>::aux_row1, start_col, N);
}

template<typename T>
inline const SubViewRow<T> SubViewRow<T>::tail(const size_t N) const
{
	const size_t start_col = SubView<T>::aux_col1 + SubView<T>::n_cols - N;

	return SubViewRow<T>(SubView<T>::m, SubView<T>::aux_row1, start_col, N);
}

template<typename T>
inline SubViewRow<T>::SubViewRow(const Matrix<T>& in_m, const size_t in_row) : SubView<T>(in_m, in_row, 0, 1, in_m.getCols())
{

}

template<typename T>
inline SubViewRow<T>::SubViewRow(const Matrix<T>& in_m, const size_t in_col, const size_t in_row1, const size_t in_n_cols) :
	SubView<T>(in_m, in_row1, in_col, 1, in_n_cols)
{

}

template<typename T>
inline std::ostream & operator<<(std::ostream & stream, const SubView<T>& mat)
{
	for (size_t i = 0;i < mat.n_rows;++i)
	{
		for (size_t f = 0;f < mat.n_cols;++f)
		{
			stream << mat(i, f) << ' ';
		}
		stream << '\n';
	}
	return stream;
}

#endif // !SUBVIEW_HPP