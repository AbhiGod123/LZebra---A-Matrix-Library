#ifndef STANDALONE_H
#define STANDALONE_H

#include "Include.h"

namespace LZebra {
	namespace random {
		double gaussianRandom();
		double uniformFloatRandom();
		int uniformIntRandom(int min = 0, int max = INT_MAX);
	}

	namespace tenseopr {
		templ Matrix<T> abs(cmat m);
		templ T accu(cmat m);
		templ Matrix<T> affmul(cmat m1, noncmat m2);

		templ Matrix<char> all(cmat m, uchar c = 0);
		templ Matrix<char> any(cmat m, uchar c = 0);
		templ bool approx_equal(cmat m1, cmat m2, uchar c, T t1, T t2 = 0);
		ftempldec Matrix<T> round(cmat m1, T tolerance);

		template<typename C, typename T> Matrix<C> arg(cpmat complexmat);
		templ T as_scalar(cmat m);
		templ Matrix<T> clamp(cmat m, T min, T max);

		templ Matrix<std::complex<T>> conj(cpmat complexmat);
		template<typename C, typename T> Matrix<C> conv_to(cmat m);
		template<typename C, typename T> Matrix<std::complex<C>> conv_to(cpmat m);

		templ Matrix<T> cross(cmat m1, cmat m2);
		templ Matrix<T> cumsum(cmat m, size_t dim = 0);
		templ Matrix<T> cumprod(cmat m, size_t dim = 0);

		templ Matrix<double> ref(cmat m, std::string ones = "ones", T* mp = nullptr);
		templ Matrix<double> rref(cmat m, std::string ones = "ones");
		templ double det(cmat m);
		templ Matrix<T> diagmat(cmat m, int val = 0);
		templ Matrix<T> diagvec(cmat m, int val = 0);

		templ Matrix<T> diff(cmat m1, size_t k = 1, uchar dim = 0);
		templ double dot(ccvec v1, ccvec v2);
		templ double norm_dot(cmat v1, cmat v2);
		templ double magnitude(cmat v1);

		templ Matrix<size_t> find(cmat m, size_t k = 0, uchar s = 0);
		templ Matrix<size_t> find_finite(cmat m);
		templ Matrix<size_t> find_nonfinite(cmat m);
		templ Matrix<size_t> find_unique(cmat m, bool ascending = 1);

		templ Matrix<T> fliplr(cmat m);
		templ Matrix<T> flipud(cmat m);
		templ Matrix<T> real(cpmat m);
		templ Matrix<T> imag(cpmat m);
		templ ColVector<size_t> ind2sub(cmat size, size_t index);
		templ Matrix<size_t> ind2sub(size_t rows, size_t cols, cmat indices);
		templ Matrix<size_t> ind2sub(cmat size, cmat indices);

		templ Matrix<size_t> index_max(cmat m, uchar dim = 0);
		templ Matrix<size_t> index_min(cmat m, uchar dim = 0);
		templ void inplace_trans(noncmat m);
		templ void inplace_trans(noncpmat m);
		templ void inplace_strans(noncpmat m);

		templ bool is_finite(cmat m);
		templ Matrix<T> join_rows(cmat m1, cmat m2);
		templ Matrix<T> join_cols(cmat m1, cmat m2);
		templ Matrix<T> join_horiz(cmat m1, cmat m2);
		templ Matrix<T> join_vert(cmat m1, cmat m2);

		templ Matrix<T> kron(cmat m1, cmat m2);

		templ Matrix<T> max(cmat m, uchar dim = 0);
		templ Matrix<T> min(cmat m, uchar dim = 0);
		templ ColVector<T> nonzeros(cmat m);

		templ Matrix<T> prod(cmat m, uchar dim = 0);
		templ size_t rank(cmat m);

		templ Matrix<T> reshape(cmat m, size_t n_rows, size_t n_cols);
		templ Matrix<T> resize(cmat m, size_t n_rows, size_t n_cols);

		templ Matrix<T> reverse(cmat m, size_t dim = 0);
		templ Matrix<T> shift(cmat m, int c, uchar dim = 0);
		templ Matrix<T> shuffle(cmat m, uchar dim = 0);
		templ ColVector<T> sort(ccvec m, std::string type = "ascend");

		templ Matrix<T> sum(cmat m, uchar dim = 0);
		templ size_t sub2ind(cmat size, size_t i, size_t j);
		templ ColVector<size_t> sub2ind(size_t rows, size_t cols, cmat indices);
		templ ColVector<size_t> sub2ind(cmat size, cmat indices);

		templ T trace(cmat m);
		templ Matrix<T> trans(cmat m);
		templ Matrix<std::complex<T>> trans(cpmat m);
		templ Matrix<std::complex<T>> strans(cpmat m);
		templ Matrix<T> trapz(cmat m, uchar dim = 0);

		templ Matrix<T> trimatu(cmat m, int k = 0);
		templ Matrix<T> trimatl(cmat m, int k = 0);
		templ Matrix<T> unique(cmat m);
		templ Matrix<T> vectorise(cmat m, uchar dim = 0);

		_fPtr Matrix<double> misc(cmat m);

		//DECOMPOSITION
		templ void lu(noncmat l, noncmat u, cmat m);
		templ Matrix<double> chol(cmat m);
		templ Matrix<double> inv(cmat m);
		templ Matrix<double> null(cmat m);

		//PROCESSING
		templ ColVector<T> conv(noncvec vec1, ccvec vec2, size_t stride = 1, std::string shape = "same");
		templ ColVector<T> conv2(noncmat mat1, cmat mat2, size_t stride = 1, std::string shape = "same");
		templ ColVector<std::complex<double>> dft(ccvec x);
		templ ColVector<std::complex<double>> dft(cpvec x);
	}
};
#endif