#ifndef STANDALONE_H
#define STANDALONE_H

#define _USE_MATH_DEFINES
#define templ template<typename T>
#define cmat const Matrix<T>& 
#define ccvec const ColVector<T>& 
#define crvec const RowVector<T>& 

template<typename T>
class Matrix;
template<typename T>
class ColVector;
template<typename T>
class RowVector;

#include <cmath>
#include <limits>

namespace random {
	float gaussianRandom();
	float uniformFloatRandom();
	int uniformIntRandom(int max = INT_MAX);
}

namespace tenseopr {
	templ Matrix<T> abs(cmat m);
	templ T accu(cmat m);
	templ T as_scalar(cmat m);
	templ Matrix<T> clamp(cmat m, T min, T max);
	templ ColVector<T> cross(ccvec v1, ccvec v2);
	templ RowVector<T> cross(crvec v1, crvec v2);
	templ Matrix<T> cumsum(cmat m, size_t dim = 0); //DOESNT WORK
	templ Matrix<T> cumprod(cmat m, size_t dim = 0); //DOESNT WORK
	templ double dot(ccvec v1, ccvec v2);
	templ double dot(crvec v1, crvec v2);
	templ double norm_dot(ccvec v1, ccvec v2);
	templ double norm_dot(crvec v1, crvec v2);
	templ double magnitude(ccvec v1);
	templ double magnitude(crvec v1);
	
}

#endif