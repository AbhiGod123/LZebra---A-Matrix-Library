#ifndef ROW_VECTOR
#define ROW_VECTOR

#include "Matrix.hpp"
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
	inline void insert(const T val);
	inline void insert_zeros(size_t s);
	inline void insert_ones(size_t s);

	//OPERATOR
	inline void operator=(const Matrix<T>& x);

	//DEFAULT CONSTRUCTOR AND DESTRUCTOR
	inline RowVector(size_t size);
	inline RowVector(size_t size, T elem);
	inline RowVector(const std::initializer_list<T>& list);
	inline RowVector(const std::vector<T>& list);
	inline RowVector(size_t size, const T* list);
	inline RowVector(const Matrix<T>& mat);
	inline RowVector();
};

#endif // !ROW_VECTOR