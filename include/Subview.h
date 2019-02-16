#ifndef SUBVIEW_H
#define SUBVIEW_H

#include "Include.h"

template<typename T>
class SubView {
public:
	const Matrix<T>& m;

	const size_t aux_row1;
	const size_t aux_col1;

	const size_t n_rows;
	const size_t n_cols;
	const size_t n_elem;

protected:
	SubView(const Matrix<T>& in_m, const size_t in_row1, const size_t in_col1, const size_t in_n_rows, const size_t in_n_cols);

public:
	inline void operator=  (const T val);
	inline void operator+= (const T val);
	inline void operator-= (const T val);
	inline void operator*= (const T val);
	inline void operator/= (const T val);

	inline void operator=  (const SubView<T>& x);
	inline void operator+= (const SubView<T>& x);
	inline void operator-= (const SubView<T>& x);
	inline void operator%= (const SubView<T>& x);

	inline void replace(const T old_val, const T new_val);

	inline void fill(const T val);
	inline void zeros();
	inline void ones();
	inline void eye();
	inline void randu();
	inline void randn();

	inline T& operator[](const size_t ii);
	inline T operator[](const size_t ii) const;

	inline T& operator()(const size_t ii);
	inline T operator()(const size_t ii) const;

	inline T& operator()(const size_t in_row, const size_t in_col);
	inline T operator()(const size_t in_row, const size_t in_col) const;

	inline T& at(const size_t in_row, const size_t in_col);
	inline T  at(const size_t in_row, const size_t in_col) const;

	inline bool is_vec()    const;
	inline bool is_finite() const;
		   
	inline bool is_inf() const;
	inline bool is_nan() const;

	inline       SubViewRow<T> row(const size_t row_num);
	inline const SubViewRow<T> row(const size_t row_num) const;

	inline       SubViewCol<T> col(const size_t col_num);
	inline const SubViewCol<T> col(const size_t col_num) const;
	
	inline       SubView<T> rows(const size_t in_row1, const size_t in_row2);
	inline const SubView<T> rows(const size_t in_row1, const size_t in_row2) const;

	inline       SubView<T> cols(const size_t in_col1, const size_t in_col2);
	inline const SubView<T> cols(const size_t in_col1, const size_t in_col2) const;

	inline       SubView<T> submat(const size_t in_row1, const size_t in_col1, const size_t in_row2, const size_t in_col2);
	inline const SubView<T> submat(const size_t in_row1, const size_t in_col1, const size_t in_row2, const size_t in_col2) const;

	private:
		friend class Matrix<T>;
		SubView();
};

template<typename T>
class SubViewCol : public SubView<T> {
public:
	inline       SubViewCol<T> subvec(const size_t in_row1, const size_t in_row2);
	inline const SubViewCol<T> subvec(const size_t in_row1, const size_t in_row2) const;

	inline       SubViewCol<T> head(const size_t N);
	inline const SubViewCol<T> head(const size_t N) const;

	inline       SubViewCol<T> tail(const size_t N);
	inline const SubViewCol<T> tail(const size_t N) const;

protected:
	friend class Matrix<T>;
	friend class ColVector<T>;
	friend class SubView<T>;

	inline SubViewCol(const Matrix<T>& in_m, const size_t in_col);
	inline SubViewCol(const Matrix<T>& in_m, const size_t in_col, const size_t in_row1, const size_t in_n_rows);
};

template<typename T>
class SubViewRow : public SubView<T> {
public:
	inline       SubViewRow<T> subvec(const size_t in_row1, const size_t in_row2);
	inline const SubViewRow<T> subvec(const size_t in_row1, const size_t in_row2) const;
						
	inline       SubViewRow<T> head(const size_t N);
	inline const SubViewRow<T> head(const size_t N) const;
						
	inline       SubViewRow<T> tail(const size_t N);
	inline const SubViewRow<T> tail(const size_t N) const;

protected:
	friend class Matrix<T>;
	friend class RowVector<T>;
	friend class SubView<T>;

	inline SubViewRow(const Matrix<T>& in_m, const size_t in_col);
	inline SubViewRow(const Matrix<T>& in_m, const size_t in_col, const size_t in_row1, const size_t in_n_rows);
};

#endif // !SUBVIEW_H