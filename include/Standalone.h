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
#define noncmat Matrix<T>&
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
	templ Matrix<T> affmul(cmat m1, noncmat m2);

	templ Matrix<char> all(cmat m, uchar c=0);
	templ Matrix<char> any(cmat m, uchar c=0);
	templ bool approx_equal(cmat m1, cmat m2, uchar c, T t1, T t2 = 0);
	ftempldec Matrix<T> round(cmat m1, T tolerance);

	template<typename C, typename T> Matrix<C> arg(cpmat complexmat);
	templ T as_scalar(cmat m);
	templ Matrix<T> clamp(cmat m, T min, T max);

	templ double cond(cmat m);//need det first
	templ Matrix<std::complex<T>> conj(cpmat complexmat);
	template<typename C, typename T> Matrix<C> conv_to(cmat m);

	templ Matrix<T> cross(cmat m1, cmat m2);
	templ Matrix<T> cumsum(cmat m, size_t dim = 0); //DOESNT WORK
	templ Matrix<T> cumprod(cmat m, size_t dim = 0); //DOESNT WORK

	templ std::tuple<T, Matrix<T>> ref(cmat m);//code later
	templ T det(cmat m);
	templ Matrix<T> diagmat(cmat m, uchar val = 0);//vector and mat
	templ Matrix<T> diagvec(cmat m, uchar val = 0);//vector and mat
 
	templ Matrix<T> diff(cmat m1, size_t k=0);//vector and mat
	templ double dot(cmat v1, cmat v2);//only vector
	templ double norm_dot(cmat v1, cmat v2);//only vector
	templ double magnitude(cmat v1);//only vector
	
	templ Matrix<T> expmat(cmat m); //later
	templ Matrix<T> expmat_sym(cmat m);//later

	/*---------------------------------------------------------------------------------------------*/
	templ Matrix<size_t> find(cmat m);
	templ Matrix<size_t> find_finite(cmat m);
	templ Matrix<size_t> find_nonfinite(cmat m);
	templ Matrix<size_t> find_unique(cmat m);

	templ Matrix<T> fliplr(cmat m);
	templ Matrix<T> flipup(cmat m);
	templ Matrix<T> real(cpmat m);
	templ Matrix<T> imag(cpmat m);



}

#endif