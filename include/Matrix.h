#ifndef MATRIX_H
#define MATRIX_H

#include "Standalone.hpp"
#include <utility>
#include <functional>
#include <iostream>
#include <fstream>
#include <vector>
#include <initializer_list>
#include <string>

template<typename T>
class Matrix
{
protected:
	size_t rows;
	size_t cols;
	size_t size;
	std::vector<T> matrix;
public:
	//GETTERS
	inline size_t getRows() const;
	inline size_t getCols() const;
	inline size_t getSize() const;

	//ELEMENT ACCESSORS
	T& operator[] (const size_t ii);
	const T& operator[] (const size_t ii) const;
	T& at(const size_t ii);
	const T& at(const size_t ii) const;
	T& operator() (const size_t ii);
	const T& operator() (const size_t ii) const;

	T& at(const size_t in_row, const size_t in_col);
	const T& at(const size_t in_row, const size_t in_col) const;
	T& operator() (const size_t in_row, const size_t in_col);
	const T& operator() (const size_t in_row, const size_t in_col) const;

	//IN-RANGE
	bool in_range(const size_t ii) const;
	bool in_range(const size_t in_row, const size_t in_col) const;

	//SAME SIZE
	bool is_same_size(const Matrix& m) const;
	bool is_same_size(const size_t in_rows, const size_t in_cols);

	//SIZE MANIPULATION
	inline void reset();
	inline void copysize(const Matrix<T>& m);

	inline void set_size(const size_t in_rows, const size_t in_cols);
	inline void resize(const size_t in_rows, const size_t in_cols);
	inline void reshape(const size_t in_rows, const size_t in_cols);

	//SET ELEMENTS
	inline void zeros();
	inline void row_zeros(const size_t in_rows);
	inline void col_zeros(const size_t in_cols);
	inline void zeros(const size_t in_rows, const size_t in_cols);

	inline void ones();
	inline void row_ones(const size_t in_rows);
	inline void col_ones(const size_t in_cols);
	inline void ones(const size_t in_rows, const size_t in_cols);

	inline void replace(const T old_val, const T new_val);
	inline void replace_row(const size_t in_rows, const T old_val, const T new_val);
	inline void replace_col(const size_t in_cols, const T old_val, const T new_val);

	inline void fill(const T val);
	inline void fill_row(const size_t in_rows, const T val);
	inline void fill_col(const size_t in_cols, const T val);
	inline void fill_rows(size_t r1, size_t r2, const T val);
	inline void fill_cols(size_t c1, size_t c2, const T val);

	//FLOATING-POINT ONLY
	inline void randu();
	inline void row_randu(const size_t in_rows);
	inline void col_randu(const size_t in_cols);
	inline void randu(const size_t in_rows, const size_t in_cols);

	inline void randn();
	inline void row_randn(const size_t in_rows);
	inline void col_randn(const size_t in_cols);
	inline void randn(const size_t in_rows, const size_t in_cols);

	//FOR INTEGERS
	inline void randi();
	inline void row_randi(const size_t in_rows);
	inline void col_randi(const size_t in_cols);
	inline void randi(const size_t in_rows, const size_t in_cols);

	typedef std::function<T(T)> lambdaT;
	typedef T(*func_p)(T);

	//FOREACH
		//LAMBDA FUNCTION
	inline void fill(const lambdaT val);
	inline void fill_row(const size_t in_rows, const lambdaT val);
	inline void fill_col(const size_t in_cols, const lambdaT val);
	inline void fill_rows(size_t r1, size_t r2, const lambdaT val);
	inline void fill_cols(size_t c1, size_t c2, const lambdaT val);
	//FUNCTOR FUNCTION
	inline void fill(const func_p val);
	inline void fill_row(const size_t in_rows, const func_p val);
	inline void fill_col(const size_t in_cols, const func_p val);
	inline void fill_rows(size_t r1, size_t r2, const func_p val);
	inline void fill_cols(size_t c1, size_t c2, const func_p val);

	//MATRIX MANIPULATION
	inline const Matrix<T> transpose();
	inline const Matrix<T> inverse();

	//OPERATIONS
		//MATRIX - MATRIX
	inline Matrix<T> operator+ (const Matrix<T>& m);
	inline Matrix<T> operator- (const Matrix<T>& m);
	inline Matrix<T> operator* (const Matrix<T>& m);
	inline Matrix<T> operator% (const Matrix<T>& m);
	inline Matrix<T> operator/ (const Matrix<T>& m);
	inline Matrix<T>& operator= (const Matrix<T>& m);
	inline Matrix<T>& operator= (Matrix<T>&& m);
	inline Matrix<T>& operator+=(const Matrix<T>& m);
	inline Matrix<T>& operator-=(const Matrix<T>& m);
	inline Matrix<T>& operator*=(const Matrix<T>& m);
	inline Matrix<T>& operator%=(const Matrix<T>& m);
	inline Matrix<T>& operator/=(const Matrix<T>& m);
	//BOOLEAN EVALUATION
	inline bool operator!=(const Matrix<T>& m) const;
	inline bool operator==(const Matrix<T>& m) const;
	inline bool operator>=(const Matrix<T>& m) const;
	inline bool operator<=(const Matrix<T>& m) const;
	inline bool operator< (const Matrix<T>& m) const;
	inline bool operator> (const Matrix<T>& m) const;

	//MATRIX - SCALAR
	inline Matrix<T>& operator= (const T val);
	inline Matrix<T>& operator+=(const T val);
	inline Matrix<T>& operator-=(const T val);
	inline Matrix<T>& operator*=(const T val);
	inline Matrix<T>& operator/=(const T val);
	inline Matrix<T>& operator++();
	inline Matrix<T>& operator--();
	inline Matrix<T> operator++(int);
	inline Matrix<T> operator--(int);

	//MATRIX CHECKS
	bool is_empty()  const;
	bool is_row_vec() const;
	bool is_col_vec() const;
	bool is_square() const;
	bool is_sorted() const;
	bool is_symmetric();
	//FLOATING-POINT ONLY
	bool is_finite() const;
	bool is_inf() const;
	bool is_nan() const;

	//MAX-MIN ELEMENTS
	inline T min() const;
	inline T max() const;
	inline size_t indexmin() const;
	inline size_t indexmax() const;

	//ITERATOR CLASSES
	typedef typename std::vector<T>::iterator              iterator;
	typedef typename std::vector<T>::const_iterator  const_iterator;

	typedef typename std::vector<T>::iterator              col_iterator;
	typedef typename std::vector<T>::const_iterator  const_col_iterator;

	class const_row_iterator;
	
	class row_iterator { //iterators through rows and single col
	public:
		inline row_iterator(const row_iterator& X);
		inline row_iterator(Matrix<T>& in_M, const size_t in_col);

		inline T& operator* ();
		inline iterator& operator& ();

		inline row_iterator& operator++();
		inline row_iterator  operator++(int);

		inline row_iterator& operator--();
		inline row_iterator  operator--(int);

		inline bool operator!=(const       row_iterator& X) const;
		inline bool operator==(const       row_iterator& X) const;
		inline bool operator!=(const const_row_iterator& X) const;
		inline bool operator==(const const_row_iterator& X) const;

		inline void print() const;

		Matrix<T>* mat;
		iterator itr;
		size_t current_col;
		size_t current_row;
	};

	class const_row_iterator {
	public:
		inline const_row_iterator(const       row_iterator& X);
		inline const_row_iterator(const const_row_iterator& X);
		inline const_row_iterator(const Matrix<T>& in_M, const size_t in_row);

		inline const T& operator*() const;

		inline const_row_iterator& operator++();
		inline const_row_iterator  operator++(int);

		inline const_row_iterator& operator--();
		inline const_row_iterator  operator--(int);

		inline bool operator!=(const       row_iterator& X) const;
		inline bool operator==(const       row_iterator& X) const;
		inline bool operator!=(const const_row_iterator& X) const;
		inline bool operator==(const const_row_iterator& X) const;

		inline void print() const;

		const Matrix<T>* mat;
		const_iterator itr;
		size_t current_row;
		size_t current_col;
	};

	//ITERATOR FUNCTIONS
	inline       iterator  begin();
	inline const_iterator  begin() const;
	inline const_iterator cbegin() const;

	inline       iterator  end();
	inline const_iterator  end() const;
	inline const_iterator cend() const;

	inline       col_iterator begin_row(const size_t row_num);
	inline const_col_iterator begin_row(const size_t row_num) const;

	inline       col_iterator end_row(const size_t row_num);
	inline const_col_iterator end_row(const size_t row_num) const;

	inline       row_iterator begin_col(const size_t col_num);
	inline const_row_iterator begin_col(const size_t col_num) const;

	inline       row_iterator end_col(const size_t col_num);
	inline const_row_iterator end_col(const size_t col_num) const;

	//SWAP
	inline void swap(Matrix<T>& m);
	inline void swap_rows(size_t r1, size_t r2);
	inline void swap_cols(size_t c1, size_t c2);

	//INSERT 
	inline void insert_rows(size_t r1, const Matrix<T>& m);
	inline void insert_cols(size_t c1, const Matrix<T>& m);

	inline void insert_row(size_t r1, const T val);
	inline void insert_col(size_t c1, const T val);

	inline void insert_rows(size_t r1, size_t r2, const T val);
	inline void insert_cols(size_t c1, size_t c2, const T val);

	inline void insert_row_zeros(size_t r1);
	inline void insert_col_zeros(size_t c1);
	inline void insert_row_ones(size_t r1);
	inline void insert_col_ones(size_t c1);

	inline void insert_rows_zeros(size_t r1, size_t r2);
	inline void insert_cols_zeros(size_t c1, size_t c2);
	inline void insert_rows_ones(size_t r1, size_t r2);
	inline void insert_cols_ones(size_t c1, size_t c2);

	//SHED 
	inline void shed_row(size_t r1);
	inline void shed_rows(size_t rbegin, size_t rend);
	inline void shed_col(size_t c1);  //DOESN'T WORK
	inline void shed_cols(size_t cbegin, size_t cend);

	//PRINT-SAVE MATRIX
	inline void print() const; 
	inline void save(const std::string& name);
	inline void load(const std::string& name);

	//DEFAULT CONSTRUCTOR AND DESTRUCTOR
	inline Matrix(size_t rows, size_t cols);
	inline Matrix(size_t rows, size_t cols, T elem); 
	inline Matrix(size_t rows, size_t cols, const std::initializer_list<T>& list);
	inline Matrix(const std::initializer_list< std::initializer_list<T> >& list); 
	inline Matrix(size_t rows, size_t cols, const std::vector<T>& list); 
	inline Matrix(const std::vector< std::vector<T> >& list); 
	inline Matrix(size_t rows, size_t cols, const T* list); 
	inline Matrix(size_t rows, size_t cols, const T * const * list);
	inline Matrix(const Matrix<T>& m);
	inline Matrix(Matrix<T>&& m); 
	inline Matrix(); 
	inline ~Matrix();
};

template<typename T>
inline std::ostream& operator<<(std::ostream& o, const Matrix<T>& m); 

template<typename T>
inline std::ostream& operator<<(std::ostream& o, const typename Matrix<T>::row_iterator& m); 

#endif // !MATRIX