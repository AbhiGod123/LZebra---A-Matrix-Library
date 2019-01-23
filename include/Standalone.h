#ifndef STANDALONE_H
#define STANDALONE_H

#define _USE_MATH_DEFINES
#define templ template<typename T>

template<typename T>
class Matrix;

#include <cmath>
#include <limits>

namespace random {
	float gaussianRandom();
	float uniformFloatRandom();
	int uniformIntRandom(int max = INT_MAX);
}

namespace tenseopr {
	templ Matrix<T> abs(const Matrix<T>& m);
	templ T accu(const Matrix<T>& m);
}

#endif