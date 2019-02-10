#ifndef SUBVIEW_H
#define SUBVIEW_H

#include "Include.h"

template<typename T>
class SubView {
	const Matrix<T>& m;

	const size_t aux_row1;
	const size_t aux_col1;

	const size_t n_rows;
	const size_t n_cols;
	const size_t n_elem;

	inline void operator=  (const T val);
	inline void operator+= (const T val);
	inline void operator-= (const T val);
	inline void operator*= (const T val);
	inline void operator/= (const T val);

	inline void operator=  (const subview& x);
	inline void operator+= (const subview& x);
	inline void operator-= (const subview& x);
	inline void operator%= (const subview& x);

	inline void replace(const T old_val, const T new_val);

	inline void fill(const T val);
	inline void zeros();
	inline void ones();
	inline void eye();
	inline void randu();
	inline void randn();

	inline T& operator[](const size_t ii);
	inline T  operator[](const size_t ii) const;

	inline T& operator()(const size_t ii);
	inline T  operator()(const size_t ii) const;

	inline T& operator()(const size_t in_row, const size_t in_col);
	inline T  operator()(const size_t in_row, const size_t in_col) const;

	inline T&         at(const size_t in_row, const size_t in_col);
	inline T          at(const size_t in_row, const size_t in_col) const;

	inline arma_warn_unused bool is_vec()    const;
	inline arma_warn_unused bool is_finite() const;

	inline arma_warn_unused bool has_inf() const;
	inline arma_warn_unused bool has_nan() const;

	inline       subview_row<T> row(const size_t row_num);
	inline const subview_row<T> row(const size_t row_num) const;

	inline       SubViewCol<T> col(const size_t col_num);
	inline const SubViewCol<T> col(const size_t col_num) const;

	inline       subview<T> rows(const size_t in_row1, const size_t in_row2);
	inline const subview<T> rows(const size_t in_row1, const size_t in_row2) const;

	inline       subview<T> cols(const size_t in_col1, const size_t in_col2);
	inline const subview<T> cols(const size_t in_col1, const size_t in_col2) const;

	inline       subview<T> submat(const size_t in_row1, const size_t in_col1, const size_t in_row2, const size_t in_col2);
	inline const subview<T> submat(const size_t in_row1, const size_t in_col1, const size_t in_row2, const size_t in_col2) const;

	inline void swap_rows(const size_t in_row1, const size_t in_row2);
	inline void swap_cols(const size_t in_col1, const size_t in_col2);

	class const_row_iterator;

	class row_iterator
	{
	public:

		inline row_iterator();
		inline row_iterator(const row_iterator& X);
		inline row_iterator(subview<T>& in_sv, const size_t in_row, const size_t in_col);

		inline arma_warn_unused T& operator* ();

		inline                  row_iterator& operator++();
		inline arma_warn_unused row_iterator  operator++(int);

		inline arma_warn_unused bool operator!=(const       row_iterator& X) const;
		inline arma_warn_unused bool operator==(const       row_iterator& X) const;
		inline arma_warn_unused bool operator!=(const const_row_iterator& X) const;
		inline arma_warn_unused bool operator==(const const_row_iterator& X) const;

		typedef std::forward_iterator_tag iterator_category;
		typedef T                        value_type;
		typedef std::ptrdiff_t            difference_type;  // TODO: not certain on this one
		typedef T*                       pointer;
		typedef T&                       reference;

		Mat<T>* M;
		T*      current_ptr;
		size_t    current_row;
		size_t    current_col;

		const size_t aux_col1;
		const size_t aux_col2_p1;
	};

	class const_row_iterator
	{
	public:

		inline const_row_iterator();
		inline const_row_iterator(const       row_iterator& X);
		inline const_row_iterator(const const_row_iterator& X);
		inline const_row_iterator(const subview<T>& in_sv, const size_t in_row, const size_t in_col);

		inline arma_warn_unused const T& operator*() const;

		inline                  const_row_iterator& operator++();
		inline arma_warn_unused const_row_iterator  operator++(int);

		inline arma_warn_unused bool operator!=(const       row_iterator& X) const;
		inline arma_warn_unused bool operator==(const       row_iterator& X) const;
		inline arma_warn_unused bool operator!=(const const_row_iterator& X) const;
		inline arma_warn_unused bool operator==(const const_row_iterator& X) const;

		typedef std::forward_iterator_tag iterator_category;
		typedef T                        value_type;
		typedef std::ptrdiff_t           difference_type;  // TODO: not certain on this one
		typedef const T*                 pointer;
		typedef const T&                 reference;

		arma_aligned const Mat<T>* M;
		arma_aligned const T*      current_ptr;
		arma_aligned       size_t    current_row;
		arma_aligned       size_t    current_col;

		arma_aligned const size_t aux_col1;
		arma_aligned const size_t aux_col2_p1;
	};

	inline       iterator  begin();
	inline const_iterator  begin() const;
	inline const_iterator cbegin() const;

	inline       iterator  end();
	inline const_iterator  end() const;
	inline const_iterator cend() const;
};

template<typename T>
class SubViewCol : public SubView<T> {
public:
	inline void operator= (const SubView<T>& x);
	inline void operator= (const SubViewCol& x);
	inline void operator= (const T val);

	inline void fill(const T val);
	inline void zeros();
	inline void ones();

	T& operator[](const size_t i);
	T  operator[](const size_t i) const;

	inline T& operator()(const size_t i);
	inline T  operator()(const size_t i) const;

	inline       SubViewCol<T> subvec(const size_t in_row1, const size_t in_row2);
	inline const SubViewCol<T> subvec(const size_t in_row1, const size_t in_row2) const;

	inline       SubViewCol<T> head(const size_t N);
	inline const SubViewCol<T> head(const size_t N) const;

	inline       SubViewCol<T> tail(const size_t N);
	inline const SubViewCol<T> tail(const size_t N) const;

	inline T min() const;
	inline T max() const;

	inline size_t index_min() const;
	inline size_t index_max() const;
};

template<typename T>
class SubViewRow : public SubView<T> {
public:
	inline void operator= (const SubView<T>& x);
	inline void operator= (const SubViewRow& x);
	inline void operator= (const T val);

	inline void fill(const T val);
	inline void zeros();
	inline void ones();

	T& operator[](const size_t i);
	T  operator[](const size_t i) const;

	inline T& operator()(const size_t i);
	inline T  operator()(const size_t i) const;

	inline       SubViewRow<T> subvec(const size_t in_row1, const size_t in_row2);
	inline const SubViewRow<T> subvec(const size_t in_row1, const size_t in_row2) const;
						
	inline       SubViewRow<T> head(const size_t N);
	inline const SubViewRow<T> head(const size_t N) const;
						
	inline       SubViewRow<T> tail(const size_t N);
	inline const SubViewRow<T> tail(const size_t N) const;

	inline T min() const;
	inline T max() const;

	inline size_t index_min() const;
	inline size_t index_max() const;
};

#endif // !SUBVIEW_H