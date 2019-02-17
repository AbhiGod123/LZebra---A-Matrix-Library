#ifndef INCLUDE_H
#define INCLUDE_H

#define _USE_MATH_DEFINES

#include <vector>
#include <math.h>
#include <cmath>
#include <complex>
#include <random>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <limits>
#include <stdlib.h>
#include <type_traits>

#define templ template<typename T>
#define cmat const Matrix<T>&
#define noncmat Matrix<T>&
#define ccvec const ColVector<T>& 
#define crvec const RowVector<T>& 
#define noncvec ColVector<T>&
#define cpmat const Matrix<std::complex<T>>&
#define noncpmat Matrix<std::complex<T>>&
#define _fPtr template<double(*_func)(double), typename T>

#define ftempldec template<typename T, typename = typename std::enable_if<std::is_floating_point<T>::value>::type>
#define itempldec template<typename T, typename = typename std::enable_if<std::is_integral<T>::value>::type>
#define templdef template<typename T, typename>

typedef unsigned char uchar;

template<typename T>
class Matrix;
template<typename T>
class ColVector;
template<typename T>
class RowVector;
template<typename T>
class SubView;

#endif // !INCLUDE_H