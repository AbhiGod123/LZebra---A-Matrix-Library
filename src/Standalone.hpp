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
	templ Matrix<T> abs(const Matrix<T>& m)
	{
		Matrix<T> mat(m);
		void(*fPtr)(int&) = [](int& val) {val = ::abs(val);};

		mat.fill(fPtr);
		return mat;
	}
	
	templ T accu(const Matrix<T>& m)
	{
		T sum = 0;
		for (typename Matrix<T>::const_iterator itr = m.begin();itr != m.end();++itr) {
			sum += *itr;
		}
		return sum;
	}


}