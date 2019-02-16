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

	inline void operator=  (const SubView& x);
	inline void operator+= (const SubView& x);
	inline void operator-= (const SubView& x);
	inline void operator%= (const SubView& x);

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

	inline void swap_rows(const size_t in_row1, const size_t in_row2);
	inline void swap_cols(const size_t in_col1, const size_t in_col2);

	//ITERATOR CLASSES
	typedef typename T*              iterator;
	typedef typename const T*  const_iterator;

	typedef typename T*              col_iterator;
	typedef typename const T*  const_col_iterator;

	class const_row_iterator;

	class row_iterator
	{
	public:
		inline row_iterator();
		inline row_iterator(const row_iterator& X);
		inline row_iterator(SubView<T>& in_sv, const size_t in_row, const size_t in_col);

		inline T& operator* ();

		inline row_iterator& operator++();
		inline row_iterator  operator++(int);

		inline bool operator!=(const       row_iterator& X) const;
		inline bool operator==(const       row_iterator& X) const;
		inline bool operator!=(const const_row_iterator& X) const;
		inline bool operator==(const const_row_iterator& X) const;

		typedef std::forward_iterator_tag iterator_category;
		typedef T                        value_type;
		typedef std::ptrdiff_t            difference_type;  // TODO: not certain on this one
		typedef T*                       pointer;
		typedef T&                       reference;

		Matrix<T>* M;
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
		inline const_row_iterator(const SubView<T>& in_sv, const size_t in_row, const size_t in_col);

		inline const T& operator*() const;

		inline const_row_iterator& operator++();
		inline const_row_iterator  operator++(int);

		inline bool operator!=(const       row_iterator& X) const;
		inline bool operator==(const       row_iterator& X) const;
		inline bool operator!=(const const_row_iterator& X) const;
		inline bool operator==(const const_row_iterator& X) const;

		typedef std::forward_iterator_tag iterator_category;
		typedef T                        value_type;
		typedef std::ptrdiff_t           difference_type;  // TODO: not certain on this one
		typedef const T*                 pointer;
		typedef const T&                 reference;

		const Matrix<T>* M;
		const T*      current_ptr;
		size_t    current_row;
		size_t    current_col;

		const size_t aux_col1;
		const size_t aux_col2_p1;
	};

	inline       iterator  begin();
	inline const_iterator  begin() const;
	inline const_iterator cbegin() const;

	inline       iterator  end();
	inline const_iterator  end() const;
	inline const_iterator cend() const;

	private:
		friend class Matrix<T>;
		SubView();
};

template<typename T>
class SubViewCol : public SubView<T> {
public:
	inline void operator= (const SubView<T>& x);

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