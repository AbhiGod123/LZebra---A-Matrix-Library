#ifndef STANDALONE_H
#define STANDALONE_H

#define _USE_MATH_DEFINES

#include <cmath>
#include <complex>
#include <algorithm>
#include <tuple>
#include <limits>
#include <type_traits>

#define templ template<typename T>
#define cmat const Matrix<T>&
#define ccvec const ColVector<T>& 
#define crvec const RowVector<T>& 
#define cpmat const Matrix<std::complex<T>>&

#define ftempldec template<typename T, typename = typename std::enable_if<std::is_floating_point<T>::value>::type>
#define ftempldef template<typename T, typename>

typedef unsigned char uchar;

template<typename T>
class Matrix;
template<typename T>
class ColVector;
template<typename T>
class RowVector;

namespace random {
	double gaussianRandom();
	double uniformFloatRandom();
	int uniformIntRandom(int max = INT_MAX);
}

namespace tenseopr {
	templ Matrix<T> abs(cmat m);
	templ T accu(cmat m);
	templ bool approx_equal(cmat m1, cmat m2, uchar c, T t1, T t2=0);
	ftempldec Matrix<T> round(cmat m1, T tolerance);

	ftempldec Matrix<T> arg(cpmat complexmat);
	templ T as_scalar(cmat m);
	templ Matrix<T> clamp(cmat m, T min, T max);

	templ double cond(cmat m);
	templ Matrix<std::complex<T>> conj(cpmat complexmat);
	template<typename T, typename C> Matrix<T> conv_to(cmat m);

	templ ColVector<T> cross(ccvec v1, ccvec v2);
	templ RowVector<T> cross(crvec v1, crvec v2);
	templ Matrix<T> cumsum(cmat m, size_t dim = 0); //DOESNT WORK
	templ Matrix<T> cumprod(cmat m, size_t dim = 0); //DOESNT WORK

	templ std::tuple<T, Matrix<T>> ref(cmat m);
	templ T det(cmat m);


	templ double dot(ccvec v1, ccvec v2);
	templ double dot(crvec v1, crvec v2);
	templ double norm_dot(ccvec v1, ccvec v2);
	templ double norm_dot(crvec v1, crvec v2);
	templ double magnitude(ccvec v1);
	templ double magnitude(crvec v1);
	
}

#endif