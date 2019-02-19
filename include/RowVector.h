#ifndef ROW_VECTOR
#define ROW_VECTOR

#include "Matrix.hpp"

namespace LZebra {
	//ONLY SUPPOSED TO BE USED FOR INTEGRAL AND FLOATING TYPES
	template<typename T>
	class RowVector : public Matrix<T>
	{
	public:
		//SIZE MANIPULATION
		inline void set_size(const size_t size);
		inline void reshape(const size_t size);
		inline void resize(const size_t size);

		//INSERT
		inline void insert(const RowVector<T>& object);
		inline void insert(size_t s, const T val);
		inline void insert(size_t index, size_t s, const T val);
		inline void insert(const T val);
		inline void insert_zeros(size_t s);
		inline void insert_ones(size_t s);

		//SUBVIEW
		inline SubViewRow<T> get_row(size_t r1, size_t r2);
		inline SubViewRow<T> tail_row(size_t r1);
		inline SubViewRow<T> head_row(size_t r2);

		//ERASE
		inline void shed(size_t i1, size_t i2);
		inline void shed_head(size_t i);
		inline void shed_tail(size_t i);

		//OPERATOR
		inline void operator=(const Matrix<T>& x);
		inline void operator=(const SubView<T>& V);

		//DEFAULT CONSTRUCTOR AND DESTRUCTOR
		inline RowVector(size_t size);
		inline RowVector(size_t size, T elem);
		inline RowVector(const std::initializer_list<T>& list);
		inline RowVector(const std::vector<T>& list);
		inline RowVector(size_t size, const T* list);
		inline RowVector(const Matrix<T>& mat);
		inline RowVector(const SubView<T>& V);
		inline RowVector();
	};
};
#endif // !ROW_VECTOR