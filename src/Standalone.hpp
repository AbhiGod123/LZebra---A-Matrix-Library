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
					return 0;
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
			T sum = 1;
			for (typename Matrix<T>::iterator itr = mat.begin(); itr != mat.end();++itr) {
				sum *= (*itr);
				*itr = sum;
			}
		}
		else if (!dim) {
			for (size_t i = 0;i < m.getRows();++i) {
				typename Matrix<T>::col_iterator itrend = mat.end_row(i);

				T cm = 1;
				for (typename Matrix<T>::col_iterator itr1 = mat.begin_row(i);itr1 != itrend;++itr1) {
					cm *= *itr1;
					*itr1 = cm;
				}
			}
		}

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
		
		return mat;
	}

	templ Matrix<T> ref(cmat m1)
	{
		Matrix<double> mat = tenseopr::conv_to<double>(m1);

		size_t m = m1.getRows();
		size_t n = m1.getCols();

		size_t pivot = 0;
		for (size_t i = 0;i < n-1;++i) {	

			bool same_col = 1;

			for (size_t c = i + pivot * m;c < mat.getSize();c+=m) {
					
				if (mat(c) != 0)
				{
					same_col = 0;
					break;
				}		
			}

			if (same_col) //check if elem range also contains the value
			{
				continue;
			}
			else {
				//we know the col doesn't contain all zeros. So there must be at least 1 number to use for the pivot point
				size_t prev_pivot = pivot;

				for (size_t j = prev_pivot;j < m;++j) {
					if (mat(j, i))
					{
						pivot = j;
						break;
					}
				}//finds the pivot point

				mat.swap_rows(prev_pivot, pivot);
				pivot = prev_pivot;

				double cur_elem = mat(prev_pivot, i);

				for (size_t row = prev_pivot + 1;row < m;++row) {
					//the goal is to make everything under the pivot zeros
					double multiplier = -1 * mat(row, i) / cur_elem;

					if (multiplier == 0)
						continue;

					for (size_t col = i;col < n;++col) {
						mat(row, col) += mat(prev_pivot, col) * multiplier;
					}
				}
				++pivot;
			}
		}
		
		return tenseopr::conv_to<T>(mat);
	}

	templ T det(cmat m)
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


		std::tuple<T,Matrix<T>> tpref = ref(m); //performs the ref and returns that matrix with a scalar. Also checks if it's a square matrix

		for (size_t i = 0;i < std::get<1>(tpref).getRows();++i) {
			std::get<0>(tpref)*=std::get<1>(tpref)(i,i);
		}

		return std::get<0>(tpref);
	}

	templ Matrix<T> diagmat(cmat m, int val)
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

	templ Matrix<T> diagvec(cmat m, int val)
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

	templ Matrix<size_t> find_unique(cmat m, bool ascending)
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

		if (ascending)
			std::sort(reali.begin(), reali.end());

		return ColVector<size_t>(reali);
	}

	templ Matrix<T> unique(cmat m)
	{
		std::vector<T> indices(m.max() + 1);
		std::vector<T> reali;

		for (size_t i = 0;i < m.getSize();++i) {
			++indices[m(i)];
		}

		for (size_t i = 1;i <= m.max();++i) {
			if (indices[i] == 1)
				reali.emplace_back(i);
		}

		return ColVector<T>(reali);
	}

	templ Matrix<T> fliplr(cmat m)
	{
		Matrix<T> mat(m);

		for (size_t i = 0;i < m.getRows();++i) {
			std::reverse(mat.begin_row(i), mat.end_row(i));
		}
		
		return mat;
	}

	templ Matrix<T> flipud(cmat m)
	{
		Matrix<T> mat(m);

		for (size_t i = 0;i < m.getCols();++i) {
			std::reverse(mat.begin_col(i), mat.end_col(i));
		}

		return mat;
	}

	templ Matrix<T> real(cpmat m)
	{
		Matrix<T> mat(m.getRows(),m.getCols());

		for (size_t i = 0;i < m.getSize();++i) {
			mat(i) = m(i).real();
		}

		return mat;
	}

	templ Matrix<T> imag(cpmat m)
	{
		Matrix<T> mat(m.getRows(), m.getCols());

		for (size_t i = 0;i < m.getSize();++i) {
			mat(i) = m(i).imag();
		}

		return mat;
	}

	templ ColVector<size_t> ind2sub(cmat size, size_t index)
	{
		ColVector<size_t> subscripts(2);

		subscripts(0) = index / size.getCols();
		subscripts(1) = index % size.getCols();

		return subscripts;
	}

	templ Matrix<size_t> ind2sub(size_t rows, size_t cols, cmat indices)
	{
		Matrix<size_t> mat(indices.getSize(), 2);

		for (size_t i = 0;i < indices.getSize();++i) {
			mat(i, 0) = indices(i) / cols;
			mat(i, 1) = indices(i) % cols;
		}

		return mat;
	}

	templ Matrix<size_t> ind2sub(cmat size, cmat indices)
	{
		Matrix<size_t> mat(indices.getSize(), 2);

		for (size_t i = 0;i < indices.getSize();++i) {
			mat(i, 0) = indices(i) / size.getCols();
			mat(i, 1) = indices(i) % size.getCols();
		}

		return mat;
	}

	templ Matrix<size_t> index_max(cmat m, uchar dim)
	{
		if (m.is_vec()) {
			return m.indexmax();
		}
		else if(!dim){
			RowVector<size_t> index(m.getCols());

			for (size_t i = 0;i < m.getCols();++i) {
				index(i) = m.indexmax_col(i);
			}

			return index;
		}
		ColVector<size_t> index(m.getRows());

		for (size_t i = 0;i < m.getRows();++i) {
			index(i) = m.indexmax_row(i);
		}

		return index;
	}

	templ Matrix<size_t> index_min(cmat m, uchar dim)
	{
		if (m.is_vec()) {
			return m.indexmin();
		}
		else if (!dim) {
			RowVector<size_t> index(m.getCols());

			for (size_t i = 0;i < m.getCols();++i) {
				index(i) = m.indexmin_col(i);
			}

			return index;
		}
		ColVector<size_t> index(m.getRows());

		for (size_t i = 0;i < m.getRows();++i) {
			index(i) = m.indexmin_row(i);
		}

		return index;
	}

	templ void inplace_trans(noncmat m)
	{
		Matrix<T> trans = m.transpose();
		m = std::move(trans);
	}

	templ void inplace_trans(noncpmat m)
	{
		Matrix<std::complex<T>> trans = m.transpose();
		m = std::move(tenseopr::conj(trans));
	}

	templ void inplace_strans(noncpmat m)
	{
		Matrix<std::complex<T>> trans = m.transpose();
		m = std::move(trans);
	}

	templ bool is_finite(cmat m)
	{
		return m.is_finite();
	}

	templ Matrix<T> join_rows(cmat m1, cmat m2)
	{
		Matrix<T> mat(m1);
		mat.insert_rows(m2);

		return mat;
	}

	templ Matrix<T> join_cols(cmat m1, cmat m2)
	{
		Matrix<T> mat(m1);
		mat.insert_cols(m2);

		return mat;
	}

	templ Matrix<T> join_horiz(cmat m1, cmat m2)
	{
		return tenseopr::join_rows(m1, m2);
	}

	templ Matrix<T> join_vert(cmat m1, cmat m2)
	{
		return tenseopr::join_cols(m1, m2);
	}

	templ Matrix<T> kron(cmat m1, cmat m2)
	{
		Matrix<T> mat(m1.getRows()*m2.getRows(), m1.getCols()*m2.getCols());

		for (size_t i = 0;i < mat.getRows();++i) {
			for (size_t j = 0;j < mat.getCols();++j) {
				size_t m1x = i / m2.getRows();
				size_t m1y = j / m2.getCols();

				size_t m2x =  i % m2.getRows();
				size_t m2y =  j % m2.getCols();
				mat(i, j) = m1(m1x, m1y) * m2(m2x, m2y);
			}
		}

		return mat;
	}

	templ Matrix<T> max(cmat m, uchar dim)
	{
		if (m.is_vec()) {
			return m.max();
		}
		else if (!dim) {
			RowVector<size_t> val(m.getCols());

			for (size_t i = 0;i < m.getCols();++i) {
				val(i) = m.max_col(i);
			}

			return val;
		}
		ColVector<size_t> val(m.getRows());

		for (size_t i = 0;i < m.getRows();++i) {
			val(i) = m.max_row(i);
		}

		return val;
	}

	templ Matrix<T> min(cmat m, uchar dim)
	{
		if (m.is_vec()) {
			return m.min();
		}
		else if (!dim) {
			RowVector<size_t> val(m.getCols());

			for (size_t i = 0;i < m.getCols();++i) {
				val(i) = m.min_col(i);
			}

			return val;
		}
		ColVector<size_t> val(m.getRows());

		for (size_t i = 0;i < m.getRows();++i) {
			val(i) = m.min_row(i);
		}

		return val;
	}

	templ ColVector<T> nonzeros(cmat m)
	{
		std::vector<T> values;
		values.reserve(m.getSize());

		for (auto itr = m.begin();itr != m.end();++itr) {
			if ((*itr))
				values.emplace_back(*itr);
		}

		return ColVector<T>(std::move(values));
	}

	templ Matrix<T> prod(cmat m, uchar dim)
	{
		if (m.is_vec()) {
			return std::accumulate(m.begin(), m.end(), 1, std::multiplies<T>());
		}
		else if (!dim) {
			ColVector<T> mat(m.getRows());

			for (size_t i = 0;i < m.getRows();++i)
				mat(i) = std::accumulate(m.begin_row(i), m.end_row(i), 1, std::multiplies<T>());

			return mat;
		}
		RowVector<T> mat(m.getRows());

		for (size_t i = 0;i < m.getRows();++i)
			mat(i) = std::accumulate(m.begin_col(i), m.end_col(i), 1, std::multiplies<T>());

		return mat;
	}

	templ size_t rank(cmat m)
	{
		size_t rank = 0;

		Matrix<T> echelon = tenseopr::ref(m);

		for (size_t i = 0;i < m.getRows();++i) {
			if (!echelon.is_same_row(i, 0))
				++rank;
		}

		return rank;
	}

	templ Matrix<T> repelem(cmat m, size_t rowcopy, size_t colcopy)
	{
		Matrix<T> mat(m);

		for (size_t i = 0;i < rowcopy;++i) {
			mat.insert_rows(mat.getRows(), mat);
		}

		for (size_t i = 0;i < colcopy;++i) {
			mat.insert_cols(mat.getCols(), mat);
		}

		return mat;
	}

	templ Matrix<T> repmat(cmat m, size_t rowcopy, size_t colcopy)
	{
		Matrix<T> mat(m);

		

		return mat;
	}

	templ Matrix<T> reshape(cmat m, size_t n_rows, size_t n_cols)
	{
		Matrix<T> mat(m);
		mat.resize(n_rows, n_cols);

		return mat;
	}

	templ Matrix<T> reverse(cmat m, size_t dim)
	{	
		if (m.is_vec()){
			Matrix<T> mat(m);
			std::reverse(mat.begin(), mat.end());

			return mat;
		}
		else if(!dim){
			Matrix<T> mat(m);

			for (size_t i = 0;i < mat.getRows();++i) {
				std::reverse(mat.begin_row(i),mat.end_row(i));
			}

			return mat;
		}
		Matrix<T> mat(m);

		for (size_t i = 0;i < mat.getCols();++i) {
			std::reverse(mat.begin_col(i), mat.end_col(i));
		}

		return mat;
	}

	templ Matrix<T> shift(cmat m, int c, uchar dim)
	{
		if (!c)
			return m;

		if (m.is_vec()) {
			Matrix<T> mat(m);

			if (c > 0)
				std::rotate(mat.begin(), mat.begin() + mat.getSize() - c, mat.end());
			else
				std::rotate(mat.begin(), mat.begin() - c, mat.end());

			return mat;
		}
		else if (!dim) {
			Matrix<T> mat(m);

			if (c > 0) {
				for (size_t i = 0;i < mat.getRows();++i) {
					std::rotate(mat.begin_row(i), mat.begin_row(i) + mat.getCols() - c, mat.end_row(i));
				}
			}
			else {
				for (size_t i = 0;i < mat.getRows();++i) {
					std::rotate(mat.begin_row(i), mat.begin_row(i) - c, mat.end_row(i));
				}
			}

			return mat;
		}

		Matrix<T> mat(m);

		if (c > 0) {
			for (size_t i = 0;i < mat.getCols();++i) {
				std::rotate(mat.begin_col(i), mat.begin_col(i) + (mat.getRows() - c), mat.end_col(i));
			}
		}
		else {
			for (size_t i = 0;i < mat.getCols();++i) {
				std::rotate(mat.begin_col(i), mat.begin_col(i) - c, mat.end_col(i));
			}
		}

		return mat;
	}

	templ Matrix<T> shuffle(cmat m, uchar dim)
	{
		std::random_device rd;
		std::mt19937 g(rd());

		if (m.is_vec()) {
			Matrix<T> mat(m);
			std::shuffle(mat.begin(), mat.end(), g);

			return mat;
		}
		else if (!dim) {
			Matrix<T> mat(m);

			for(size_t i=0;i<mat.getRows();++i)
				std::shuffle(mat.begin_row(i), mat.end_row(i), g);
			
			return mat;
		}
		Matrix<T> mat(m);

		for (size_t i = 0;i < mat.getCols();++i)
			std::shuffle(mat.begin_col(i), mat.end_col(i), g);
		
		return mat;
	}

	templ Matrix<T> sort(cmat m, std::string type, uchar dim)
	{
		if (m.is_vec()) {
			Matrix<T> mat(m);
			if(type == "ascend")
			std::sort(mat.begin(), mat.end());
			else
				std::sort(mat.begin(), mat.end(), std::greater<T>());

			return mat;
		}
		else if (!dim) {
			Matrix<T> mat(m);

			for (size_t i = 0;i < mat.getCols();++i){
				if (type == "ascend")
					std::sort(mat.begin_row(i), mat.end_row(i));
				else
					std::sort(mat.begin_row(i), mat.end_row(i), std::greater<T>());
			}

			return mat;
		}
		Matrix<T> mat(m);

		for (size_t i = 0;i < mat.getRows();++i) {
			if (type == "ascend")
				std::sort(mat.begin_col(i), mat.end_col(i));
			else
				std::sort(mat.begin_col(i), mat.end_col(i), std::greater<T>());
		}

		return mat;
	}

	templ Matrix<T> sum(cmat m, uchar dim)
	{
		if (m.is_vec()) {
			return std::accumulate(m.begin(), m.end(), 0);
		}
		else if (!dim) {
			ColVector<T> mat(m.getRows());

			for (size_t i = 0;i < m.getRows();++i)
				mat(i) = std::accumulate(m.begin_row(i), m.end_row(i), 0);

			return mat;
		}
		RowVector<T> mat(m.getRows());

		for (size_t i = 0;i < m.getRows();++i)
			mat(i) = std::accumulate(m.begin_col(i), m.end_col(i), 0);

		return mat;
	}

	templ size_t sub2ind(cmat size, size_t i, size_t j)
	{
		return (i * size.getCols() + j);
	}

	templ ColVector<size_t> sub2ind(size_t rows, size_t cols, cmat indices)
	{
		ColVector<size_t> vec(indices.getRows());

		for (size_t i = 0;i < indices.getRows();++i) {
			vec(i) = indices(i, 0) * cols + indices(i, 1);
		}

		return vec;
	}

	templ ColVector<size_t> sub2ind(cmat size, cmat indices)
	{
		ColVector<size_t> vec(indices.getRows());

		for (size_t i = 0;i < indices.getRows();++i) {
			vec(i) = indices(i, 0) * size.getCols() + indices(i, 1);
		}

		return vec;
	}

	templ T trace(cmat m)
	{
		T sum = 0;
		for (size_t i = 0;i < m.getCols() && i < m.getRows();++i) {
			sum += mat(i, i);
		}

		return sum;
	}

	templ Matrix<T> trans(cmat m)
	{
		Matrix<T> mat(m.getCols(), m.getRows());

		for (size_t i = 0;i < m.getCols();++i) {
			for (size_t f = 0;f < m.getRows();++f) {
				mat(i, f) = m(f, i);
			}
		}

		return mat;
	}

	templ Matrix<std::complex<T>> trans(cpmat m)
	{
		Matrix<std::complex<T>> trans(m.getCols(), m.getRows());

		for (size_t i = 0;i < m.getCols();++i) {
			for (size_t f = 0;f < m.getRows();++f) {
				trans(i, f) = m(f, i);
			}
		}

		return tenseopr::conj(trans);
	}

	templ Matrix<std::complex<T>> strans(cpmat m)
	{
		Matrix<std::complex<T>> trans(m.getCols(), m.getRows());

		for (size_t i = 0;i < m.getCols();++i) {
			for (size_t f = 0;f < m.getRows();++f) {
				trans(i, f) = m(f, i);
			}
		}

		return trans;
	}

	templ Matrix<T> trapz(cmat m, uchar dim)
	{
		if (!dim) {
			ColVector<T> mat(m.getRows());

			for (size_t i = 0;i < m.getRows();++i) {
				T sum = m(i,0);

				for (size_t f = 1;f < m.getCols()-1;++f) {
					sum += 2 * m(i, f);
				}
				sum += m(i, m.getCols() - 2);

				mat(i) = 0.5 * sum;
			}

			return mat;
		}
		RowVector<T> mat(m.getCols());

		for (size_t i = 0;i < m.getCols();++i) {
			T sum = m(0, i);

			for (size_t f = 1;f < m.getRows() - 1;++f) {
				sum += 2 * m(f, i);
			}
			sum += m(m.getCols() - 2,i);

			mat(i) = 0.5 * sum;
		}

		return mat;
	}

	templ Matrix<T> trimatu(cmat m, int val)
	{
		Matrix<T> mat;
		mat.copysize(m);

		if (val == 0) {
			for (size_t i = 0;i < m.getRows()-m.getCols()-2;++i) {
				for (size_t f = 1;f < mat.getCols();++f) {
					mat(i, f) = m(i, f);
				}
			}
			return mat;
		}

		if (val > 0) {
			for (size_t i = 0;i < m.getCols() - val - 1;++i) {
				for (size_t f = val + i + 1;f < mat.getCols();++f) {
					mat(i, f) = m(i, f);
				}
			}

			return mat;
		}

		for (size_t i = 0;i < m.getCols() - val && i < m.getRows();++i) {
			for (size_t f = i;f <m.getCols();++f) {
				mat(i, f) = m(i, f);
			}
		}

		return mat;
	}

	templ Matrix<T> trimatl(cmat m, int val)
	{
		Matrix<T> mat;
		mat.copysize(m);

		if (val == 0) {
			for (size_t i = 1;i < m.getRows();++i) {
				for (size_t f = 0;f < val + i && f < mat.getCols();++f) {
					mat(i, f) = m(i, f);
				}
			}
			return mat;
		}

		if (val > 0) {
			for (size_t i = 0;i < m.getRows();++i) {
				for (size_t f = 0;f < val + i && f < mat.getCols();++f) {
					mat(i, f) = m(i, f);
				}
			}
			
			return mat;
		}

		for (size_t i = -val + 1;i < m.getRows();++i) {
			for (size_t f = 0;f < 1 + i - (-val + 1) && f < mat.getCols();++f) {
				mat(i, f) = m(i, f);
			}
		}

		return mat;
	}

	templ Matrix<T> vectorise(cmat m, uchar dim)
	{
		if (!dim) {
			return tenseopr::reshape(m, m.getSize(), 1);
		}
		return tenseopr::reshape(m, 1, m.getSize());
	}

	_fPtr Matrix<double> misc(cmat m)
	{
		Matrix<double> mat = tenseopr::conv_to<double>(m);

		mat.fill(_func);
		return mat;
	}

}