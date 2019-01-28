#include "Standalone.h"

namespace random {
	double gaussianRandom() {
		double v1, v2, s;
		do {
			v1 = 2.0f * static_cast<double>(rand() % 10000) / (10000) - 1.0f;
			v2 = 2.0f * static_cast<double>(rand() % 10000) / (10000) - 1.0f;
			s = v1 * v1 + v2 * v2;
		} while (s >= 1.0f || s == 0.0f);

		s = pow((-2.0f * log(s)) / s, 0.5f);

		return v1 * s;
	}

	double uniformFloatRandom()
	{
		return static_cast<double>(rand() % 10000) / (10000);
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

	templ Matrix<T> affmul(cmat m1, noncmat m2)
	{
		m2.insert_row(m2.getRows(), 1);

		Matrix<T> resultmat = const_cast<Matrix<T>&>(m1) * m2;
		m2.shed_row(m2.getRows()-1);

		return resultmat;
	}

	templ bool all(cmat m)
	{
		for (auto itr = m.begin(); itr != m.end();++itr) {
			if (!(*itr))
				return false;
		}
		return true;
	}

	templ bool any(cmat m)
	{
		for (auto itr = m.begin(); itr != m.end();++itr) {
			if (*itr)
				return true;
		}
		return false;
	}

	templ bool approx_equal(cmat m1, cmat m2, uchar c, T t1, T t2)
	{
		if (m1.getCols() != m2.getCols() || m1.getRows() != m2.getRows())
			return false;

		switch (c)
		{
			case 0:
			{
				for (size_t i = 0;i < m1.getSize();++i) {
					if (::abs(m1(i) - m2(i)) > t1)
						return false;
				}
			}
			break;
			case 1:
			{
				for (size_t i = 0;i < m1.getSize();++i) {
					if (::abs(m1(i) - m2(i))/std::max(::abs(m1(i)), ::abs(m2(i))) > t1)
						return false;
				}
			}
			break;
			case 2:
			{
				for (size_t i = 0;i < m1.getSize();++i)
				if (::abs(m1(i) - m2(i)) > t1 && ::abs(m1(i) - m2(i)) / std::max(::abs(m1(i)), ::abs(m2(i))) > t2)
					return false;
			}
			break;
		}
		return true;
		
	}

	ftempldef Matrix<T> round(cmat m1, T tol)
	{
		Matrix<T> appmat(m1);

		for (auto itr = appmat.begin();itr != appmat.end();++itr) {
			if (::abs(::round(*itr) - *itr) < tol)
				*itr = ::round(*itr);
		}

		return appmat;
	}

	ftempldef Matrix<T> arg(cpmat complexmat)
	{
		Matrix<T> phaseangles;
		phaseangles.copysize(complexmat);

		for (size_t i = 0;i < phaseangles.getSize();++i) {
			phaseangles(i) = ::atan2(complexmat.imag(), complexmat.real());
		}
		
		return phaseangles;
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

		uchar x = 0;
		uchar y = 1;
		uchar z = 2;

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

	templ std::tuple<T, Matrix<T>> ref(cmat m)
	{
		std::tuple<T, Matrix<T>> badtuple = std::make_tuple(0, Matrix<T>());

		if (m.getCols() != m.getRows())
		{
			std::cout << "Non-square matrix!" << std::endl;
			return std::make_tuple(0, Matrix<T>());
		}

		for (size_t i = 0;i < m.getRows();++i)
			if (m.is_same_row(i, 0) || m.is_same_col(i, 0))
				return badtuple;
		//A. IF AN ENTIRE ROW OR COL CONTAINS ZERO THEN DET = 0

		for (size_t i = 0;i < m.getRows() - 1;++i) 
			for (size_t k = i + 1;k < m.getRows();++k) 
				if (m.is_equal_rows(i, k) || m.is_equal_cols(i, k))
					return badtuple;
		
		//B. IF ANY 2 ROWS OR COLS ARE EQUAL THEN DET = 0
		
		
		//C. IF 2 ROWS OR COLS ARE PROPORTIONAL TO EACH OTHER THEN DET = 0



		T value = 1;
		Matrix<T> reducedmat(m);
		
				

		return std::make_tuple(value, reducedmat);
	}

	templ T det(cmat m)
	{
		std::tuple<T,Matrix<T>> tpref = ref(m); //performs the rref and returns that matrix with a scalar. Also checks if it's a square matrix

		for (size_t i = 0;i < std::get<1>(tpref).getRows();++i) {
			std::get<0>(tpref)*=std::get<1>(tpref)(i,i);
		}

		return std::get<0>(tpref);
	}

	templ Matrix<T> diagmat(cmat m, uchar val)
	{
		Matrix<T> dmat;
		dmat.copysize(m);

		if (val >= dmat.getSize())
		{
			std::cout << "Value over size" << std::endl;
		}

		for (size_t i = 0;i < dmat.getSize();++i) {
			dmat(i, i+val) = m(i, i+val);
		}
		return dmat;
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