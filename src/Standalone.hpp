#include "Standalone.h"

namespace random {
	float gaussianRandom() {
		float v1, v2, s;
		do {
			v1 = 2.0f * (float)(rand() % 10000) / (10000) - 1.0f;
			v2 = 2.0f * (float)(rand() % 10000) / (10000) - 1.0f;
			s = v1 * v1 + v2 * v2;
		} while (s >= 1.0f || s == 0.0f);

		s = (float)pow((-2.0f * log(s)) / s, 0.5f);

		return v1 * s;
	}

	float uniformFloatRandom()
	{
		return static_cast<float>(rand() % 10000) / (10000);
	}

	int uniformIntRandom(int max)
	{
		return rand() % max;
	}
}

namespace tenseopr
{
	templ Matrix<T> abs(cmat m)
	{
		Matrix<T> mat(m);
		void(*fPtr)(int&) = [](int& val) {val = ::abs(val);};

		mat.fill(fPtr);
		return mat;
	}
	
	templ T accu(cmat m)
	{
		T sum = 0;
		for (typename Matrix<T>::const_iterator itr = m.begin();itr != m.end();++itr) {
			sum += *itr;
		}
		return sum;
	}

	templ T as_scalar(cmat m) {
		return m[0];
	}

	templ Matrix<T> clamp(cmat m,T min, T max) {
		Matrix<T> mat(m);

		for (typename Matrix<T>::iterator itr = mat.begin();itr != mat.end();++itr) {
			if (*itr < min)
				*itr = min;
			else if (*itr > max)
				*itr = max;
		}
		return mat;
	}

	templ ColVector<T> cross(ccvec v1, ccvec v2) {
		ColVector<T> vec(3);

		size_t x = 0;
		size_t y = 1;
		size_t z = 2;

		vec(x) = v1(y)*v2(z) - v1(z)*v2(y);
		vec(y) = v1(z)*v2(x) - v1(x)*v2(z);
		vec(z) = v1(x)*v2(y) - v1(y)*v2(x);

		return vec;
	}

	templ RowVector<T> cross(crvec v1, crvec v2) {
		RowVector<T> vec(3);

		size_t x = 0;
		size_t y = 1;
		size_t z = 2;

		vec(x) = v1(y)*v2(z) - v1(z)*v2(y);
		vec(y) = v1(z)*v2(x) - v1(x)*v2(z);
		vec(z) = v1(x)*v2(y) - v1(y)*v2(x);

		return vec;
	}

	templ Matrix<T> cumsum(cmat m, size_t dim)
	{
		Matrix<T> mat(m);

		switch (dim) {
			case 0: {
				for (size_t i = 0;i < m.getRows();++i) {
					typename Matrix<T>::col_iterator itrend = mat.end_row(i);

					T cm = 0;
					for (typename Matrix<T>::col_iterator itr1 = mat.begin_row(i);itr1 != itrend;++itr1) {
						cm += *itr1;
						*itr1 = cm;
					}
				}
				break;
			}
			case 1: {
				for (size_t i = 0;i < m.getCols();++i) {
					typename Matrix<T>::row_iterator itrend = mat.end_col(i);
					
					T cm = 0;
					for (typename Matrix<T>::row_iterator itr1 = mat.begin_col(i); itr1 != itrend;++itr1) {
						cm += *itr1;
						std::cout << cm << std::endl;
						*itr1 = cm;
						
					}
				}
				break;
			}
		}

		return mat;
	}

	templ Matrix<T> cumprod(cmat m, size_t dim)
	{
		Matrix<T> mat;
		mat.copysize(m);

		switch (dim) {
			case 0: {
				for (size_t i = 0;i < m.getRows();++i) {
					typename Matrix<T>::col_iterator itrend = mat.end_row(i);

					T cm = 0;
					for (typename Matrix<T>::col_iterator itr1 = mat.begin_row(i);itr1 != itrend;++itr1) {
						cm *= *itr1;
						*itr1 = cm;
					}
				}
				break;
			}
			case 1: {
				for (size_t i = 0;i < m.getCols();++i) {
					typename Matrix<T>::row_iterator itrend = mat.end_col(i);

					T cm = 0;
					for (typename Matrix<T>::row_iterator itr1 = mat.begin_col(i); itr1 != itrend;++itr1) {
						cm *= *itr1;
						std::cout << cm << std::endl;
						*itr1 = cm;

					}
				}
				break;
			}
		}

		return mat;
	}

	templ double dot(ccvec v1, ccvec v2)
	{
		double sum = 0;
		for (size_t i = 0;i < v1.getSize();++i) {
			sum+= v1(i) * v2(i);
		}
		return sum;
	}

	templ double dot(crvec v1, crvec v2) 
	{
		double sum = 0;
		for (size_t i = 0;i < v1.getSize();++i) {
			sum += v1(i) * v2(i);
		}
		return sum;
	}

	templ double norm_dot(ccvec v1, ccvec v2)
	{
		return dot(v1, v2) / (magnitude(v1) * magnitude(v2));
	}

	templ double norm_dot(crvec v1, crvec v2)
	{
		return dot(v1, v2) / (magnitude(v1) * magnitude(v2));
	}

	templ double magnitude(ccvec v1)
	{
		double sum = 0;

		for (size_t i = 0;i < v1.getSize();++i) {
			sum += ::pow(v1(i),2);
		}
		return pow(sum,0.5);
	}

	templ double magnitude(crvec v1)
	{
		double sum = 0;

		for (size_t i = 0;i < v1.getSize();++i) {
			sum += ::pow(v1(i), 2);
		}
		return pow(sum, 0.5);
	}

	
}