#ifndef COL_VECTOR
#define COL_VECTOR

#include "Matrix.hpp"
template<typename T>
class ColVector : public Matrix<T>
{
public:
	//SIZE MANIPULATION
	inline void set_size(const size_t size);
	inline void reshape(const size_t size);
	inline void resize(const size_t size);

	//INSERT
	inline void insert(const ColVector<T>& object);
	inline void insert(size_t s, const T val);
	inline void insert(const T val);
	inline void insert_zeros(size_t s);
	inline void insert_ones(size_t s);

	//DEFAULT CONSTRUCTOR AND DESTRUCTOR
	inline ColVector(size_t size);
	inline ColVector(size_t size, T elem);
	inline ColVector(const std::initializer_list<T>& list);
	inline ColVector(const std::vector<T>& list);
	inline ColVector(size_t size, const T* list);
	inline ColVector();
};

#endif // !COL_VECTOR