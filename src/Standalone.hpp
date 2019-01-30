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

	templ Matrix<char> all(cmat m, uchar c)
	{
		if (m.is_vec()) {
			for (auto itr = m.begin(); itr != m.end();++itr) {
				if (!(*itr))
					return Matrix<char>(1,1,0);
			}
		}
		else if (!c) {
			RowVector<char> c(m.getCols());

			for (size_t i = 0;i < m.getCols();++i) {

				typename Matrix<T>::row_iterator itrend = m->end_col(i);

				for (typename Matrix<T>::row_iterator itr = m->begin_col(i); itr != itrend;++itr) {
					if (!(*itr)) {
						c(i) = 0;
						break;
					}
				}
			}

			return c;
		}
		else {
			ColVector<char> c(m.getRows());

			for (size_t i = 0;i < m.getRows();++i) {

				typename Matrix<T>::col_iterator itrend = m->end_row(i);

				for (typename Matrix<T>::col_iterator itr = m->begin_row(i); itr != itrend;++itr) {
					if (!(*itr)) {
						c(i) = 0;
						break;
					}
				}
			}

			return c;
		}
		
		return Matrix<char>(1, 1, 1);
	}

	templ Matrix<char> any(cmat m, uchar c)
	{
		if (m.is_vec()) {
			for (auto itr = m.begin(); itr != m.end();++itr) {
				if ((*itr))
					return Matrix<char>(1, 1, 1);
			}
		}
		else if (!c) {
			RowVector<char> c(m.getCols());

			for (size_t i = 0;i < m.getCols();++i) {

				typename Matrix<T>::row_iterator itrend = m->end_col(i);

				for (typename Matrix<T>::row_iterator itr = m->begin_col(i); itr != itrend;++itr) {
					if ((*itr)) {
						c(i) = 1;
						break;
					}
				}
			}

			return c;
		}
		else {
			ColVector<char> c(m.getRows());

			for (size_t i = 0;i < m.getRows();++i) {

				typename Matrix<T>::col_iterator itrend = m->end_row(i);

				for (typename Matrix<T>::col_iterator itr = m->begin_row(i); itr != itrend;++itr) {
					if ((*itr)) {
						c(i) = 1;
						break;
					}
				}
			}

			return c;
		}

		return Matrix<char>(1, 1, 0);
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

	template<typename C, typename T>
	Matrix<C> arg(cpmat complexmat)
	{
		Matrix<C> phaseangles(complexmat.getRows(),complexmat.getCols());

		for (size_t i = 0;i < phaseangles.getSize();++i) {
			phaseangles(i) = ::atan2(complexmat(i).imag(), complexmat(i).real());
		}

		return phaseangles;
	}

	templ T as_scalar(cmat m) {
		return m[0];
	}

	templ Matrix<T> clamp(cmat m,T min, T max) {
		Matrix<T> mat(m);

		for (auto itr = mat.begin();itr != mat.end();++itr) {
			if (*itr < min)
				*itr = min;
			else if (*itr > max)
				*itr = max;
		}
		return mat;
	}

	templ Matrix<std::complex<T>> conj(cpmat complexmat)
	{
		Matrix<std::complex<T>> conjugate(complexmat);

		for (auto itr = conjugate.begin();itr != conjugate.end();++itr) {
			(*itr).imag(-1 * (*itr).imag());
		}
		return conjugate;
 	}

	template<typename C, typename T>
	Matrix<C> conv_to(cmat m)
	{
		Matrix<C> casted(m.getRows(),m.getCols());

		for (size_t i = 0;i < casted.getSize();++i) {
			casted(i) = static_cast<C>(m(i));
		}
		return casted;
	}

	templ Matrix<T> cross(cmat m1, cmat m2)
	{
		Matrix<T> mat(3);

		constexpr char x = 0;
		constexpr char y = 1;
		constexpr char z = 2;

		mat(x) = m1(y)*m2(z) - m1(z)*m2(y);
		mat(y) = m1(z)*m2(x) - m1(x)*m2(z);
		mat(z) = m1(x)*m2(y) - m1(y)*m2(x);

		return mat;
	}

	templ Matrix<T> cumsum(cmat m, size_t dim)
	{
		Matrix<T> mat(m);

		if (m.is_vec()) {
			T sum = 0;
			for (typename Matrix<T>::iterator itr = mat.begin(); itr != mat.end();++itr) {
				sum += (*itr);
				*itr = sum;
			}
		}
		else {
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
					T early = mat(0);

					for (size_t i = 0;i < m.getCols();++i) {
						typename Matrix<T>::row_iterator itrend = mat.end_col(i);

						T cm = 0;
						for (typename Matrix<T>::row_iterator itr1 = mat.begin_col(i); itr1 != itrend;++itr1) {
							cm += *itr1;
							*itr1 = cm;
						}
					}
					mat(0) = early;
					break;
				}
			}
		}
		
		return mat;
	}

	templ Matrix<T> cumprod(cmat m, size_t dim)
	{
		Matrix<T> mat(m);

		if (m.is_vec()) {
			T sum = 0;
			for (typename Matrix<T>::iterator itr = mat.begin(); itr != mat.end();++itr) {
				sum *= (*itr);
				*itr = sum;
			}
		}
		else {
			switch (dim) {
				case 0: {
					for (size_t i = 0;i < m.getRows();++i) {
						typename Matrix<T>::col_iterator itrend = mat.end_row(i);

						T cm = 1;
						for (typename Matrix<T>::col_iterator itr1 = mat.begin_row(i);itr1 != itrend;++itr1) {
							cm *= *itr1;
							*itr1 = cm;
						}
					}
					break;
				}
				case 1: {
					T early = mat(0);

					for (size_t i = 0;i < m.getCols();++i) {
						typename Matrix<T>::row_iterator itrend = mat.end_col(i);
						
						T cm = 1;
						for (typename Matrix<T>::row_iterator itr1 = mat.begin_col(i); itr1 != itrend;++itr1) {
							cm *= *itr1;
							*itr1 = cm;
						}
					}
					
					mat(0) = early;
					break;
				}
			}
		}

		return mat;
	}

	templ std::tuple<T, Matrix<T>> ref(cmat m)
	{
		std::tuple<T, Matrix<T>> badtuple = std::make_tuple(0, Matrix<T>());

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
		std::tuple<T,Matrix<T>> tpref = ref(m); //performs the ref and returns that matrix with a scalar. Also checks if it's a square matrix

		for (size_t i = 0;i < std::get<1>(tpref).getRows();++i) {
			std::get<0>(tpref)*=std::get<1>(tpref)(i,i);
		}

		return std::get<0>(tpref);
	}

	templ Matrix<T> diagmat(cmat m, char val)
	{
		Matrix<T> dmat;

		if (m.is_vec()) {
			dmat.set_size(m.getSize() + ::abs(val), m.getSize() + ::abs(val));
			std::cout << dmat.getSize() << std::endl;

			if (val>=0) {
				for (size_t i = 0;i < m.getSize();++i) {
					dmat(i, i + val) = m(i);
				}
			}
			else {
				for (size_t i = 0;i < m.getSize();++i) {
					dmat(i-val, i) = m(i);
				}
			}
		}
		else {
			dmat.copysize(m);

			if (val >= 0) {
				for (size_t i = 0;i < m.getCols() - val;++i) {
					dmat(i, i+val) = m(i, i+val);
				}
			}
			else {
				for (size_t i = 0;i < m.getRows() + val;++i) {
					dmat(i-val, i) = m(i-val, i);
				}
			}
		}

		return dmat;
	}

	templ Matrix<T> diagvec(cmat m, char val)
	{
		ColVector<T> vec;

		if (val >= 0) {
			vec.set_size(m.getCols() - val);

			for (size_t i = 0;i < m.getCols() - val;++i) {
				vec(i) = m(i, i + val);
			}
		}
		else {
			vec.set_size(m.getRows() + val);

			for (size_t i = 0;i < m.getRows() + val;++i) {
				vec(i) = m(i-val, i);
			}
		}
		return vec;
	}

	templ Matrix<T> diff(cmat m1, size_t k, uchar dim)
	{
		Matrix<T> mat;

		if (m1.is_vec()) {
			if (m1.getCols() == 1)
			{
				mat.set_size(1, m1.getRows() - 1);
			}
			else {
				mat.set_size(m1.getCols() - 1,1);
			}

			for (size_t i = 1;i < m1.getSize();++i) {
				mat(i - 1) = m1(i) - m1(i - 1);
			}
		}
		else {
			if (m1.getRows() == 1 || m1.getCols() == 1) {
				return mat;
			}
			if (!dim) {
				mat.set_size(m1.getRows() - 1, m1.getCols());

				for (size_t i = 0;i < mat.getCols();++i) {
					for (size_t f = 1;f < m1.getRows();++f) {
						mat(f - 1,i) = m1(f,i) - m1(f - 1,i);
					}
				}
			}
			else {
				mat.set_size(m1.getRows(), m1.getCols() - 1);

				for (size_t i = 0;i < mat.getRows();++i) {
					for (size_t f = 1;f < m1.getCols();++f) {
						mat(i, f - 1) = m1(i, f) - m1(i, f - 1);
					}
				}
			}
		}

		if (k > 1)
			return tenseopr::diff(mat, k - 1, dim);
		return mat;
	}

	templ double dot(cmat v1, cmat v2)
	{
		double sum = 0;

		if (!v1.is_vec() || !v2.is_vec())
		{
			std::cout << "Not a vector" << '\n';
		}
		else {
			for (size_t i = 0;i < v1.getSize();++i) {
				sum += v1(i) * v2(i);
			}
		}

		return sum;
	}

	templ double norm_dot(cmat v1, cmat v2)
	{
		return dot(v1, v2) / (magnitude(v1) * magnitude(v2));
	}

	templ double magnitude(cmat v1)
	{
		double sum = 0;

		if (!v1.is_vec())
		{
			std::cout << "Not a vector" << '\n';
		}
		else {
			for (size_t i = 0;i < v1.getSize();++i) {
				sum += ::pow(v1(i), 2);
			}
		}

		return pow(sum,0.5);
	}

	templ Matrix<size_t> find(cmat m, size_t k, uchar s)
	{
		std::vector<size_t> indices;
		indices.reserve(m.getSize() / 2);
		
		if (!k) { //RETURN ALL NON-ZERO INDICES  
			for (size_t i = 0;i < m.getSize();++i) {
				if(m(i))
				indices.emplace_back(i);
			}
		}
		else {
			if (!s) {
				for (size_t i = 0;i < m.getSize();++i) {
					if (m(i))
						indices.emplace_back(i);

					if (indices.size() == k)
						break;
				}
			}
			else {
				for (size_t i = m.getSize()-1;i >= 0;--i) {
					if (m(i))
						indices.emplace_back(i);

					if (indices.size() == k)
						break;
				}
			}
		}

		return ColVector<size_t>(std::move(indices));
	}

	templ Matrix<size_t> find_finite(cmat m)
	{
		std::vector<size_t> indices;
		indices.reserve(m.getSize() / 2);

		for (size_t i = 0;i < m.getSize();++i) {
			if (m.is_finite())
				indices.emplace_back(i);
		}

		return ColVector<size_t>(std::move(indices));
	}

	templ Matrix<size_t> find_nonfinite(cmat m)
	{
		std::vector<T> indices;
		indices.reserve(m.getSize() / 2);

		for (size_t i = 0;i < m.getSize();++i) {
			if (!m.is_finite())
				indices.emplace_back(i);
		}

		return ColVector<size_t>(std::move(indices));
	}

	templ Matrix<size_t> find_unique(cmat m)
	{
		std::vector<size_t> indices(m.max() + 1);
		std::vector<size_t> reali;

		for (size_t i = 0;i < m.getSize();++i) {
			++indices[m(i)];
		}

		for (size_t i = 1;i <= m.max();++i) {
			if(indices[i]==1)
			reali.emplace_back(i);
		}

		return ColVector<size_t>(reali);
	}


}