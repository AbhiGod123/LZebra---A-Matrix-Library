#include "Matrix.h"

template<typename T>
inline T & Matrix<T>::operator[](size_t index)
{
	return matrix[index];
}

template<typename T>
const T & Matrix<T>::operator[](const size_t ii) const
{
	return matrix[ii];
}

template<typename T>
T & Matrix<T>::at(const size_t ii)
{
	return matrix[ii];
}

template<typename T>
const T & Matrix<T>::at(const size_t ii) const
{
	return matrix[ii];
}

template<typename T>
T & Matrix<T>::operator()(const size_t ii)
{
	return matrix[ii];
}

template<typename T>
const T & Matrix<T>::operator()(const size_t ii) const
{
	return matrix[ii];
}

template<typename T>
T & Matrix<T>::at(const size_t in_row, const size_t in_col)
{
	return matrix[in_row * cols + in_col];
}

template<typename T>
const T & Matrix<T>::at(const size_t in_row, const size_t in_col) const
{
	return matrix[in_row * cols + in_col];
}

template<typename T>
T & Matrix<T>::operator()(const size_t in_row, const size_t in_col)
{
	return matrix[in_row * cols + in_col];
}

template<typename T>
const T & Matrix<T>::operator()(const size_t in_row, const size_t in_col) const
{
	return matrix[in_row * cols + in_col];
}

template<typename T>
bool Matrix<T>::in_range(const size_t ii) const
{
	return ii < size;
}

template<typename T>
bool Matrix<T>::in_range(const size_t in_row, const size_t in_col) const
{
	return this->in_range(in_row * cols + in_col);
}

template<typename T>
inline bool Matrix<T>::is_same_size(const Matrix & m) const
{
	return rows == m.rows && cols == m.cols;
}

template<typename T>
inline bool Matrix<T>::is_same_size(const size_t in_rows, const size_t in_cols)
{
	return rows == in_rows && cols == in_cols;
}

template<typename T>
inline void Matrix<T>::reset()
{
	matrix.clear();
	rows = 0;
	cols = 0;
	size = 0;
}

template<typename T>
inline void Matrix<T>::copysize(const Matrix<T>& m)
{
	this->set_size(m.rows, m.cols);
}

template<typename T>
inline void Matrix<T>::set_size(const size_t in_rows, const size_t in_cols)
{
	this->reset();
	matrix = std::vector<T>();
	rows = in_rows;
	cols = in_cols;
	size = in_rows * in_cols;
	matrix.resize(size);
}

template<typename T>
inline void Matrix<T>::resize(const size_t in_rows, const size_t in_cols)
{
	rows = in_rows;
	cols = in_cols;
	size = in_rows * in_cols;
	matrix.resize(size);
}

template<typename T>
inline void Matrix<T>::zeros()
{
	this->fill(0);
}

template<typename T>
inline void Matrix<T>::row_zeros(const size_t in_rows)
{
	this->fill_row(in_rows, 0);
}

template<typename T>
inline void Matrix<T>::col_zeros(const size_t in_cols)
{
	this->fill_col(in_cols, 0);
}

template<typename T>
inline void Matrix<T>::zeros(const size_t in_rows, const size_t in_cols)
{
	this->set_size(in_rows, in_cols);
	for (T& i : matrix) {
		i = 0;
	}
}

template<typename T>
inline void Matrix<T>::ones()
{
	this->fill(1);
}

template<typename T>
inline void Matrix<T>::row_ones(const size_t in_rows)
{
	this->fill_row(in_rows, 1);
}

template<typename T>
inline void Matrix<T>::col_ones(const size_t in_cols)
{
	this->fill_col(in_cols, 1);
}

template<typename T>
inline void Matrix<T>::ones(const size_t in_rows, const size_t in_cols)
{
	this->set_size(in_rows, in_cols);
	for (T& i : matrix) {
		i = 1;
	}
}

template<typename T>
inline void Matrix<T>::replace(const T old_val, const T new_val)
{
	for (T& i : matrix) {
		if (i == old_val)
			i = new_val;
	}
}

template<typename T>
inline void Matrix<T>::replace_row(const size_t in_rows, const T old_val, const T new_val)
{
	Matrix<T>::col_iterator itrend = end_row(in_rows);

	for (Matrix<T>::col_iterator itr = begin_row(in_rows); itr != itrend;itr++) {
		if (*itr == old_val)
			*itr = new_val;
	}
}

template<typename T>
inline void Matrix<T>::replace_col(const size_t in_cols, const T old_val, const T new_val)
{
	Matrix<T>::row_iterator itrend = end_col(in_cols);

	for (Matrix<T>::row_iterator itr = begin_col(in_cols); itr != itrend;itr++) {
		if (*itr == old_val)
			*itr = new_val;
	}
}

template<typename T>
inline void Matrix<T>::fill(const T val)
{
	for (T& i : matrix) {
		i = val;
	}
}

template<typename T>
inline void Matrix<T>::fill_row(const size_t in_rows, const T val)
{
	Matrix<T>::col_iterator itr = begin_row(in_rows);
	Matrix<T>::col_iterator itrend = end_row(in_rows);

	std::fill(itr, itrend, val);
}

template<typename T>
inline void Matrix<T>::fill_col(const size_t in_cols, const T val)
{
	Matrix<T>::row_iterator itrend = end_col(in_cols);

	for (Matrix<T>::row_iterator itr = begin_col(in_cols); itr != itrend;itr++) {
		*itr = val;
	}
}

template<typename T>
inline void Matrix<T>::fill_rows(size_t r1, size_t r2, const T val)
{
	for (size_t i = r1;i <= r2;++i) {
		this->fill_row(i, val);
	}
}

template<typename T>
inline void Matrix<T>::fill_cols(size_t c1, size_t c2, const T val)
{
	for (size_t i = c1;i <= c2;++i) {
		this->fill_col(i, val);
	}
}

template<typename T>
inline void Matrix<T>::randu()
{
	for (T& i : matrix) {
		i = uniformFloatRandom();
	}
}

template<typename T>
inline void Matrix<T>::row_randu(const size_t in_rows)
{
	typename std::vector<T>::iterator begin = matrix.begin() + in_rows * cols;
	typename std::vector<T>::iterator end = begin + cols;

	for (auto it = begin; it != end; it++) {
		*it = uniformFloatRandom();
	}
}

template<typename T>
inline void Matrix<T>::col_randu(const size_t in_cols)
{
	typename std::vector<T>::iterator begin = matrix.begin() + in_cols;
	typename std::vector<T>::iterator end = begin * rows;

	for (auto it = begin; it != end; it += rows) {
		*it = uniformFloatRandom();
	}
}

template<typename T>
inline void Matrix<T>::randu(const size_t in_rows, const size_t in_cols)
{
	this->set_size(in_rows, in_cols);
	for (T& i : matrix) {
		 i = uniformFloatRandom();
	}
}

template<typename T>
inline void Matrix<T>::randn()
{
	for (T& i : matrix) {
		i = gaussianRandom();
	}
}

template<typename T>
inline void Matrix<T>::row_randn(const size_t in_rows)
{
	typename std::vector<T>::iterator begin = matrix.begin() + in_rows * cols;
	typename std::vector<T>::iterator end = begin + cols;

	for (auto it = begin; it != end; it++) {
		*it = gaussianRandom();
	}
}

template<typename T>
inline void Matrix<T>::col_randn(const size_t in_cols)
{
	typename std::vector<T>::iterator begin = matrix.begin() + in_cols;
	typename std::vector<T>::iterator end = begin * rows;

	for (auto it = begin; it != end; it += rows) {
		*it = gaussianRandom();
	}
}

template<typename T>
inline void Matrix<T>::randn(const size_t in_rows, const size_t in_cols)
{
	this->set_size(in_rows, in_cols);
	for (auto i : matrix) {
		matrix[i] = gaussianRandom();
	}
}

template<typename T>
inline void Matrix<T>::randi()
{
	for (T& i : matrix) {
		i = uniformIntRandom();
	}
}

template<typename T>
inline void Matrix<T>::row_randi(const size_t in_rows)
{
	typename std::vector<T>::iterator begin = matrix.begin() + in_rows * cols;
	typename std::vector<T>::iterator end = begin + cols;

	for (auto it = begin; it != end; it++) {
		*it = uniformIntRandom();
	}
}

template<typename T>
inline void Matrix<T>::col_randi(const size_t in_cols)
{
	typename std::vector<T>::iterator begin = matrix.begin() + in_cols;
	typename std::vector<T>::iterator end = begin * rows;

	for (auto it = begin; it != end; it += rows) {
		*it = uniformIntRandom();
	}
}

template<typename T>
inline void Matrix<T>::randi(const size_t in_rows, const size_t in_cols)
{
	this->set_size(in_rows, in_cols);
	for (T& i : matrix) {
		i = uniformIntRandom();
	}
}

template<typename T>
inline void Matrix<T>::fill(const lambdaT val)
{
	for (T& i : matrix) {
		i = val(i);
	}
}

template<typename T>
inline void Matrix<T>::fill_row(const size_t in_rows, const lambdaT val)
{
	typename std::vector<T>::iterator begin = matrix.begin() + in_rows * cols;
	typename std::vector<T>::iterator end = begin + cols;

	for (auto it = begin; it != end; it++) {
		*it = val(*it);
	}
}

template<typename T>
inline void Matrix<T>::fill_col(const size_t in_cols, const lambdaT val)
{
	typename std::vector<T>::iterator begin = matrix.begin() + in_cols;
	typename std::vector<T>::iterator end = begin * rows;

	for (auto it = begin; it != end; it += rows) {
		*it = val(*it);
	}
}

template<typename T>
inline void Matrix<T>::fill_rows(size_t r1, size_t r2, const lambdaT val)
{
	for (size_t i = r1;i <= r2;++i) {
		this->fill_row(i, val);
	}
}

template<typename T>
inline void Matrix<T>::fill_cols(size_t c1, size_t c2, const lambdaT val)
{
	for (size_t i = c1;i <= c2;++i) {
		this->fill_col(i, val);
	}
}

template<typename T>
inline void Matrix<T>::fill(const func_p val)
{
	for (T& i : matrix) {
		i = val(i);
	}
}

template<typename T>
inline void Matrix<T>::fill_row(const size_t in_rows, const func_p val)
{
	typename std::vector<T>::iterator begin = matrix.begin() + in_rows * cols;
	typename std::vector<T>::iterator end = begin + cols;

	for (auto it = begin; it != end; it++) {
		*it = val(*it);
	}
}

template<typename T>
inline void Matrix<T>::fill_col(const size_t in_cols, const func_p val)
{
	typename std::vector<T>::iterator begin = matrix.begin() + in_cols;
	typename std::vector<T>::iterator end = begin * rows;

	for (auto it = begin; it != end; it += rows) {
		*it = val(*it);
	}
}

template<typename T>
inline void Matrix<T>::fill_rows(size_t r1, size_t r2, const func_p val)
{
	for (size_t i = r1;i <= r2;++i) {
		this->fill_row(i, val);
	}
}

template<typename T>
inline void Matrix<T>::fill_cols(size_t c1, size_t c2, const func_p val)
{
	for (size_t i = c1;i <= c2;++i) {
		this->fill_col(i, val);
	}
}

template<typename T>
inline const Matrix<T> Matrix<T>::transpose()
{
	Matrix<T> newmat(cols, rows);

	for (size_t i = 0;i < cols;++i) {
		for (size_t f = 0;f < rows;++f) {
			//newmat(i, f) = matrix[f * cols + i];
			newmat(i, f) = (*this)(f,i);
		}
	}

	return newmat;
}

template<typename T>
inline const Matrix<T> Matrix<T>::inverse()
{
	// TODO: insert return statement here
}

template<typename T>
inline Matrix<T> Matrix<T>::operator+(const Matrix<T>& m)
{
	Matrix<T> result(*this);
	result += m;

	return result;
}

template<typename T>
inline Matrix<T> Matrix<T>::operator-(const Matrix<T>& m)
{
	Matrix<T> result(*this);
	result -= m;

	return result;
}

template<typename T>
inline Matrix<T> Matrix<T>::operator*(const Matrix<T>& m)
{
	Matrix<T> newmat(rows, m.cols);

	for (size_t i = 0;i < rows;++i) {
		for (size_t f = 0;f < m.cols;++f) {
			for (size_t k = 0;k < cols;++k) {
				newmat(i, f) += (*this)[i * cols + k] * m(k, f);
			}
		}
	}

	return newmat;
}

template<typename T>
inline Matrix<T> Matrix<T>::operator%(const Matrix<T>& m)
{
	Matrix<T> result(*this);
	result %= m;

	return result;
}

template<typename T>
inline Matrix<T> Matrix<T>::operator/(const Matrix<T>& m)
{
	Matrix<T> result(*this);
	result /= m;

	return result;
}

template<typename T>
inline Matrix<T>& Matrix<T>::operator=(const Matrix<T>& m)
{
	set_size(m.rows, m.cols);

	for (size_t i = 0;i < m.size;i++) {
		matrix[i] = m[i];
	}

	return *this;
}

template<typename T>
inline Matrix<T>& Matrix<T>::operator=(Matrix<T>&& m)
{
	matrix = std::move(m.matrix);
	rows = std::exchange(m.rows, 0);
	cols = std::exchange(m.cols, 0);
	size = std::exchange(m.size, 0);
}

template<typename T>
inline Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& m)
{
	for (size_t i = 0;i < matrix.size();++i) {
		matrix[i] += m[i];
	}

	return *this;
}

template<typename T>
inline Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& m)
{
	for (size_t i = 0;i < matrix.size();++i) {
		matrix[i] -= m[i];
	}

	return *this;
}

template<typename T>
inline Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& m)
{
	Matrix<T> result = (*this) * m;
	(*this) = result;
	return *this;
}

template<typename T>
inline Matrix<T>& Matrix<T>::operator%=(const Matrix<T>& m)
{
	for (size_t i = 0;i < matrix.size();++i) {
		matrix[i] *= m[i];
	}

	return *this;
}

template<typename T>
inline Matrix<T>& Matrix<T>::operator/=(const Matrix<T>& m)
{
	for (size_t i = 0;i < matrix.size();++i) {
		matrix[i] /= m[i];
	}

	return *this;
}

template<typename T>
inline bool Matrix<T>::operator!=(const Matrix<T>& m) const
{
	if (!is_same_size(m))
		return false;

	for (size_t i = 0;i < matrix.size();++i) {
		if (matrix[i] == m[i])
			return false;
	}
	return true;
}

template<typename T>
inline bool Matrix<T>::operator==(const Matrix<T>& m) const
{
	if (!is_same_size(m))
		return false;

	for (size_t i = 0;i < matrix.size();++i) {
		if (matrix[i] != m[i])
			return false;
	}
	return true;
}

template<typename T>
inline bool Matrix<T>::operator>=(const Matrix<T>& m) const
{
	if (!is_same_size(m))
		return false;

	for (size_t i = 0;i < matrix.size();++i) {
		if (matrix[i] < m[i])
			return false;
	}

	return true;
}

template<typename T>
inline bool Matrix<T>::operator<=(const Matrix<T>& m) const
{
	if (!is_same_size(m))
		return false;

	for (size_t i = 0;i < matrix.size();++i) {
		if (matrix[i] > m[i])
			return false;
	}

	return true;
}

template<typename T>
inline bool Matrix<T>::operator<(const Matrix<T>& m) const
{
	if (!is_same_size(m))
		return false;

	for (size_t i = 0;i < matrix.size();++i) {
		if (matrix[i] >= m[i])
			return false;
	}

	return true;
}

template<typename T>
inline bool Matrix<T>::operator>(const Matrix<T>& m) const
{
	if (!is_same_size(m))
		return false;

	for (size_t i = 0;i < matrix.size();++i) {
		if (matrix[i] <= m[i])
			return false;
	}

	return true;
}

template<typename T>
inline Matrix<T>& Matrix<T>::operator=(const T val)
{
	set_size(1, 1);
	fill(val);

	return *this;
}

template<typename T>
inline Matrix<T>& Matrix<T>::operator+=(const T val)
{
	for (T& i : matrix) {
		i += val;
	}

	return *this;
}

template<typename T>
inline Matrix<T>& Matrix<T>::operator-=(const T val)
{
	for (T& i : matrix) {
		i -= val;
	}

	return *this;
}

template<typename T>
inline Matrix<T>& Matrix<T>::operator*=(const T val)
{
	for (T& i : matrix) {
		i *= val;
	}

	return *this;
}

template<typename T>
inline Matrix<T>& Matrix<T>::operator/=(const T val)
{
	for (T& i : matrix) {
		i /= val;
	}

	return *this;
}

template<typename T>
inline Matrix<T>& Matrix<T>::operator++()
{
	(*this) += 1;

	return *this;
}

template<typename T>
inline Matrix<T>& Matrix<T>::operator--()
{
	(*this) -= 1;

	return *this;
}

template<typename T>
bool Matrix<T>::is_empty() const
{
	return static_cast<bool>(size);
}

template<typename T>
bool Matrix<T>::is_row_vec() const
{
	return rows == 1;
}

template<typename T>
bool Matrix<T>::is_col_vec() const
{
	return cols == 1;
}

template<typename T>
bool Matrix<T>::is_square() const
{
	return rows == cols;
}

template<typename T>
bool Matrix<T>::is_sorted() const
{
	for (auto itr = matrix.begin() + 1; itr != matrix.end(); ++itr) {
		if (*(itr-1) > *itr)
			return false;
	}

	return true;
}

template<typename T>
bool Matrix<T>::is_symmetric()
{
	return (*this) == transpose();
}

template<typename T>
bool Matrix<T>::is_finite() const
{
	for (T& i : matrix) {
		if (!isfinite(i))
			return false;
	}
	return true;
}

template<typename T>
bool Matrix<T>::is_inf() const
{
	for (T& i : matrix) {
		if (!isinf(i))
			return false;
	}
	return true;
}

template<typename T>
bool Matrix<T>::is_nan() const
{
	for (T& i : matrix) {
		if (!isnan(i))
			return false;
	}
	return true;
}

template<typename T>
inline T Matrix<T>::min() const
{
	return matrix[this->indexmin()];
}

template<typename T>
inline T Matrix<T>::max() const
{
	return matrix[this->indexmax()];
}

template<typename T>
inline size_t Matrix<T>::indexmin() const
{
	T minval = std::numeric_limits<T>::max();
	size_t index = 0;

	for (size_t i = 0;i < size;++i) {
		if (matrix[i] < minval)
		{
			index = i;
			minval = matrix[i];
		}
	}

	return index;
}

template<typename T>
inline size_t Matrix<T>::indexmax() const
{
	T maxval = std::numeric_limits<T>::min();

	size_t index = 0;

	for (size_t i = 0;i < size;++i) {
		if (matrix[i] > maxval)
		{
			index = i;
			maxval = matrix[i];
		}
	}

	return index;
}

template<typename T>
inline Matrix<T>::row_iterator::row_iterator() :mat(NULL), itr(NULL), current_col(0), current_row(0)
{

}
	
template<typename T>
inline Matrix<T>::row_iterator::row_iterator(const row_iterator& X) : mat(X.mat), itr(X.itr), current_col(X.current_col), current_row(X.current_row)
{

}

template<typename T>
inline Matrix<T>::row_iterator::row_iterator(Matrix<T>& in_M, const size_t in_col) : mat(&in_M), itr(in_M.begin() + in_col), current_col(in_col), current_row(0)
{

}

template<typename T>
inline T & Matrix<T>::row_iterator::operator*()
{
	return (*itr);
}

template<typename T>
inline typename Matrix<T>::row_iterator& Matrix<T>::row_iterator::operator++()
{
	current_row++;

	if (current_row == mat->getRows())
	{
		current_row = 0;
		current_col++;

		if (current_col == mat->getCols())
			current_col = 0;

		itr = mat->begin() + current_col;
	}
	else
	{
		itr += mat->getCols();
	}

	return *this;
}

template<typename T>
inline typename Matrix<T>::row_iterator Matrix<T>::row_iterator::operator++(int)
{
	typename Matrix<T>::row_iterator temp(*this);

	++(*this);

	return temp;
}

template<typename T>
inline typename Matrix<T>::row_iterator& Matrix<T>::row_iterator::operator--()
{
	current_row--;

	if (current_row == -1) {
		current_row = mat->getRows() - 1;

		current_col--;

		if (current_col == -1)
			current_col = mat->getCols() - 1;

		itr = mat->begin() + (current_row * mat->getCols() + current_col);
	}
	else {
		itr -= mat->getCols();
	}

	return *this;
}

template<typename T>
inline typename Matrix<T>::row_iterator Matrix<T>::row_iterator::operator--(int)
{
	typename Matrix<T>::row_iterator temp(*this);

	--(*this);

	return temp;
}

template<typename T>
inline bool Matrix<T>::row_iterator::operator!=(const row_iterator & X) const
{
	return itr != X.itr;
}

template<typename T>
inline bool Matrix<T>::row_iterator::operator==(const row_iterator & X) const
{
	return itr == X.itr;
}

template<typename T>
inline bool Matrix<T>::row_iterator::operator!=(const const_row_iterator & X) const
{
	return itr != X.itr;
}

template<typename T>
inline bool Matrix<T>::row_iterator::operator==(const const_row_iterator & X) const
{
	return itr == X.itr;
}

template<typename T>
inline void Matrix<T>::row_iterator::print() const
{
	std::cout << "itr: " << (*itr) << ' ';
	std::cout << "row: " << current_row << ' ';
	std::cout << "col: " << current_col << std::endl;
}

template<typename T>
inline Matrix<T>::const_row_iterator::const_row_iterator() :mat(NULL), itr(NULL), current_col(0), current_row(0)
{

}

template<typename T>
inline Matrix<T>::const_row_iterator::const_row_iterator(const row_iterator& X) : mat(X.mat), itr(X.itr), current_col(X.current_col), current_row(X.current_row)
{

}

template<typename T>
inline Matrix<T>::const_row_iterator::const_row_iterator(const const_row_iterator& X) : mat(X.mat), itr(X.itr), current_col(X.current_col), current_row(X.current_row)
{

}

template<typename T>
inline Matrix<T>::const_row_iterator::const_row_iterator(const Matrix<T>& in_M, const size_t in_col) : mat(&in_M), itr(in_M.begin() + in_col), current_col(in_col), current_row(0)
{

}

template<typename T>
inline const T & Matrix<T>::const_row_iterator::operator*() const
{
	return (*itr);
}

template<typename T>
inline typename Matrix<T>::const_row_iterator& Matrix<T>::const_row_iterator::operator++()
{
	current_row++;

	if (current_row == mat->getRows())
	{
		current_row = 0;
		current_col++;

		if (current_col == mat->getCols())
			current_col = 0;

		itr = mat->begin() + current_col;
	}
	else
	{
		itr += mat->getCols();
	}

	return *this;
}

template<typename T>
inline typename Matrix<T>::const_row_iterator Matrix<T>::const_row_iterator::operator++(int)
{
	typename Matrix<T>::const_row_iterator temp(*this);

	++(*this);

	return temp;
}

template<typename T>
inline typename Matrix<T>::const_row_iterator& Matrix<T>::const_row_iterator::operator--()
{
	current_row--;

	if (current_row == -1) {
		current_row = mat->getRows() - 1;

		current_col--;

		if (current_col == -1)
			current_col = mat->getCols() - 1;

		itr = mat->begin() + (current_row * mat->getCols() + current_col);
	}
	else {
		itr -= mat->getCols();
	}

	return *this;
}

template<typename T>
inline typename Matrix<T>::const_row_iterator Matrix<T>::const_row_iterator::operator--(int)
{
	typename Matrix<T>::const_row_iterator temp(*this);

	--(*this);

	return temp;
}

template<typename T>
inline bool Matrix<T>::const_row_iterator::operator!=(const row_iterator & X) const
{
	return itr != X.itr;
}

template<typename T>
inline bool Matrix<T>::const_row_iterator::operator==(const row_iterator & X) const
{
	return itr == X.itr;
}

template<typename T>
inline bool Matrix<T>::const_row_iterator::operator!=(const const_row_iterator & X) const
{
	return itr != X.itr;
}

template<typename T>
inline bool Matrix<T>::const_row_iterator::operator==(const const_row_iterator & X) const
{
	return itr == X.itr;
}

template<typename T>
inline void Matrix<T>::const_row_iterator::print() const
{
	std::cout << "itr: " << (*itr) << ' ';
	std::cout << "row: " << current_row << ' ';
	std::cout << "col: " << current_col << std::endl;
}

template<typename T>
inline typename Matrix<T>::iterator Matrix<T>::begin() {
	return matrix.begin();
}

template<typename T>
inline typename Matrix<T>::const_iterator Matrix<T>::begin() const
{
	return matrix.begin();
}

template<typename T>
inline typename Matrix<T>::const_iterator Matrix<T>::cbegin() const
{
	return matrix.cbegin();
}

template<typename T>
inline typename Matrix<T>::iterator Matrix<T>::end() {
	return matrix.end();
}

template<typename T>
inline typename Matrix<T>::const_iterator Matrix<T>::end() const
{
	return matrix.end();
}

template<typename T>
inline typename Matrix<T>::const_iterator Matrix<T>::cend() const
{
	return matrix.cend();
}

template<typename T>
inline typename Matrix<T>::col_iterator Matrix<T>::begin_row(const size_t row_num)
{
	return matrix.begin() + (row_num * cols);
}

template<typename T>
inline typename Matrix<T>::const_col_iterator Matrix<T>::begin_row(const size_t row_num) const
{
	return matrix.begin() + (row_num * cols);
}

template<typename T>
inline typename Matrix<T>::col_iterator Matrix<T>::end_row(const size_t row_num)
{
	return matrix.begin() + (row_num * cols + cols);
}

template<typename T>
inline typename Matrix<T>::const_col_iterator Matrix<T>::end_row(const size_t row_num) const
{
	return matrix.begin() + (row_num * cols + cols);
}

template<typename T>
inline typename Matrix<T>::row_iterator Matrix<T>::begin_col(const size_t col_num)
{
	return row_iterator(*this,col_num);
}

template<typename T>
inline typename Matrix<T>::const_row_iterator Matrix<T>::begin_col(const size_t col_num) const
{
	return const_row_iterator(*this, col_num);
}

template<typename T>
inline typename Matrix<T>::row_iterator Matrix<T>::end_col(const size_t col_num)
{
	return row_iterator(*this, col_num + 1);
}

template<typename T>
inline typename Matrix<T>::const_row_iterator Matrix<T>::end_col(const size_t col_num) const
{
	return const_row_iterator(*this, col_num + 1);
}

template<typename T>
inline void Matrix<T>::swap(Matrix<T>& m)
{
	matrix.swap(m.matrix);
	std::swap(rows, m.rows);
	std::swap(cols, m.cols);
	std::swap(size, m.size);
}

template<typename T>
inline void Matrix<T>::swap_rows(size_t r1, size_t r2)
{
	if (r1 > r2)
		std::swap(r1, r2);
	typename std::vector<T>::iterator begin = matrix.begin() + r1 * cols;
	typename std::vector<T>::iterator end = begin + cols;

	for (auto it = begin; it != end; ++it) {
		std::swap(*it, matrix[cols * r2 + it - begin]);
	}
}

template<typename T>
inline void Matrix<T>::swap_cols(size_t c1, size_t c2)
{
	if (c1 > c2)
		std::swap(c1, c2);
	typename std::vector<T>::iterator begin = matrix.begin() + c1;
	typename std::vector<T>::iterator end = begin * rows;

	for (auto it = begin; it != end; it += rows) {
		std::swap(*it, matrix[rows * c2 + it - begin]);
	}
}

template<typename T>
inline void Matrix<T>::insert_rows(size_t r1, const Matrix<T>& m)
{
	matrix.insert(r1*cols, m.begin(), m.end());
	rows = rows + m.rows;
	size = rows * cols;
}

template<typename T>
inline void Matrix<T>::insert_cols(size_t c1, const Matrix<T>& m)
{
	typename std::vector<T>::iterator begin = matrix.begin() + c1;
	typename std::vector<T>::iterator end = begin * rows;

	for (size_t a = 0;a < m.cols;a++) {
		for (size_t i = c1;i < c1 * rows;i += rows) {
			matrix.insert(i, m(i, a));
		}
	}
	cols = cols + m.cols;
	size = rows * cols;
}

template<typename T>
inline void Matrix<T>::insert_row(size_t r1, const T val)
{
	matrix.insert(r1*cols, val, cols);
	rows++;
	size = rows * cols;
}

template<typename T>
inline void Matrix<T>::insert_col(size_t c1, const T val)
{
	typename std::vector<T>::iterator begin = matrix.begin() + c1;
	typename std::vector<T>::iterator end = begin * rows;

	for (size_t i = c1;i < c1 * rows;i += rows) {
		matrix.insert(i, val);
	}

	cols++;
	size = rows * cols;
}

template<typename T>
inline void Matrix<T>::insert_rows(size_t r1, size_t r2, const T val)
{
	for (size_t i = r1;i <= r2;++i) {
		this->insert_row(i, val);
	}
}

template<typename T>
inline void Matrix<T>::insert_cols(size_t c1, size_t c2, const T val)
{
	for (size_t i = c1;i <= c2;++i) {
		this->insert_col(i, val);
	}
}

template<typename T>
inline void Matrix<T>::insert_row_zeros(size_t r1)
{
	this->insert_row(r1, 0);
}

template<typename T>
inline void Matrix<T>::insert_col_zeros(size_t c1)
{
	this->insert_col(c1, 0);
}

template<typename T>
inline void Matrix<T>::insert_row_ones(size_t r1)
{
	this->insert_row(r1, 1);
}

template<typename T>
inline void Matrix<T>::insert_col_ones(size_t c1)
{
	this->insert_col(c1, 1);
}

template<typename T>
inline void Matrix<T>::insert_rows_zeros(size_t r1, size_t r2)
{
	this->insert_rows(r1, r2, 0);
}

template<typename T>
inline void Matrix<T>::insert_cols_zeros(size_t c1, size_t c2)
{
	this->insert_cols(c1, c2, 0);
}

template<typename T>
inline void Matrix<T>::insert_rows_ones(size_t r1, size_t r2)
{
	this->insert_rows(r1, r2, 1);
}

template<typename T>
inline void Matrix<T>::insert_cols_ones(size_t c1, size_t c2)
{
	this->insert_cols(c1, c2, 1);
}

template<typename T>
inline void Matrix<T>::shed_row(size_t r1)
{
	matrix.erase(r1*cols, r1*cols + cols);
	rows++;
	size = rows * cols;
}

template<typename T>
inline void Matrix<T>::shed_rows(size_t rbegin, size_t rend)
{
	for (size_t i = rbegin;i <= rend;++i) {
		this->shed_row(i);
	}
}

template<typename T>
inline void Matrix<T>::shed_col(size_t c1)
{
	typename std::vector<T>::iterator begin = matrix.begin() + c1;
	typename std::vector<T>::iterator end = begin * rows;


	for (size_t i = c1;i < c1 * rows;i += rows) {
		matrix.erase(i);
	}

	cols++;
	size = rows * cols;
}

template<typename T>
inline void Matrix<T>::shed_cols(size_t cbegin, size_t cend)
{
	for (size_t i = cbegin;i <= cend;++i) {
		this->shed_col(i);
	}
}

template<typename T>
inline std::ostream & operator<<(std::ostream & stream, const Matrix<T>& mat)
{
	for (size_t i = 0;i < mat.getRows();++i)
	{
		for (size_t f = 0;f < mat.getCols();++f)
		{
			stream << mat(i, f) << ' ';
		}
		stream << std::endl;
	}
	return stream;
}

template<typename T>
inline std::ostream & operator<<(std::ostream & stream, const typename Matrix<T>::row_iterator & m)
{
	stream << "itr: " << (*m.itr) << ' ';
	stream << "row: " << m.current_row << ' ';
	stream << "col: " << m.current_col << std::endl;

	return stream;
}



template<typename T>
inline void Matrix<T>::print() const
{
	for (size_t i = 0;i < rows;++i)
	{
		for (size_t f = 0;f < cols;++f)
		{
			std::cout << (*this)(i,f) << ' ';
		}
		std::cout << std::endl;
	}
}

template<typename T>
inline void Matrix<T>::save(const std::string & name)
{
	std::ofstream stream;

	stream.open(name);

	stream << rows << '\n';
	stream << cols << '\n';

	for (auto i : matrix) {
		stream << i << std::endl;
	}

	stream << std::endl;
	stream.close();
}

template<typename T>
inline void Matrix<T>::load(const std::string & name)
{
	
}

template<typename T>
inline size_t Matrix<T>::getRows() const
{
	return rows;
}

template<typename T>
inline size_t Matrix<T>::getCols() const
{
	return cols;
}

template<typename T>
inline size_t Matrix<T>::getSize() const
{
	return size;
}

template<typename T>
inline Matrix<T>& Matrix<T>::getMatrix(Matrix<T>& m)
{
	return m.matrix;
}

template<typename T>
inline Matrix<T>::Matrix(size_t rows, size_t cols) : rows(rows), cols(cols), size(rows * cols)
{
	matrix.resize(size);
}

template<typename T>
inline Matrix<T>::Matrix(size_t rows, size_t cols, T elem) : Matrix(rows, cols)
{
	for (T &i : matrix) {
		i = elem;
	}
}

template<typename T>
inline Matrix<T>::Matrix(size_t rows, size_t cols, const std::initializer_list<T>& list) : Matrix(rows, cols)
{
	matrix = list;
}

template<typename T>
inline Matrix<T>::Matrix(const std::initializer_list<std::initializer_list<T>>& list)
{
	for (auto& i : list) {
		matrix.insert(matrix.end(), i.begin(), i.end());
	}
	rows = list.size();
	cols = list.begin()[0].size();
	size = rows * cols;
}

template<typename T>
inline Matrix<T>::Matrix(size_t rows, size_t cols, const std::vector<T>& list) : Matrix(rows, cols)
{
	matrix = list;
}

template<typename T>
inline Matrix<T>::Matrix(const std::vector<std::vector<T>>& list) 
{
	for (auto& i : list) {
		matrix.insert(matrix.end(), i.begin(), i.end());
	}	
	rows = list.size();
	cols = list[0].size();
	size = rows * cols;
}

template<typename T>
inline Matrix<T>::Matrix(size_t rows, size_t cols, const T * list) : Matrix(rows, cols)
{
	for (size_t i = 0;i < size;++i) {
		matrix[i] = list[i];
	}
}

template<typename T>
inline Matrix<T>::Matrix(size_t rows, size_t cols, const T * const * list)
{
	for (size_t i = 0;i < rows;i++) {
		matrix.insert(matrix.end(), list[i], list[i] + cols);
	}
	
	this->rows = rows;
	this->cols = cols;
	this->size = rows * cols;
}

template<typename T>
inline Matrix<T>::Matrix(const Matrix<T>& m) : matrix(m.matrix), rows(m.rows), cols(m.cols), size(m.size) {
	
}

template<typename T>
inline Matrix<T>::Matrix(Matrix<T>&& m) : matrix(std::move(m.matrix)), 
rows(std::exchange(m.rows, 0)), 
cols(std::exchange(m.cols, 0)), 
size(std::exchange(m.size, 0))
{

}

template<typename T>
inline Matrix<T>::Matrix() {

}

template<typename T>
inline Matrix<T>::~Matrix() {

}