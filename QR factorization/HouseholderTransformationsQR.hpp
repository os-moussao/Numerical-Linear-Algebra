#pragma once

#include "../Matrix/Matrix.hpp"

/**
 * Getting the QR factorization of A ∈ ℝ𝑛x𝑛 using Householder Transformations
 * R = H_n-1 H_n-2 ... H_1 A , where H_i is orthogonal (transpose(H_i) = inverse(H_i)), and inverse(H_i) = H_i
 * Q = inverse(H_n-1 H_n-2 ... H_1) = H_1 H_2 ... H_n-1
 */

#define sign(x) (x < 0? -1: 1)

template<typename T>
pair<Matrix<long double>, Matrix<long double>> QR_HouseHolder(const Matrix <T> &A) {
	assert(A.rows() && A.rows()==A.cols());

	int n = A.rows();
	Matrix<long double> Q(identity(n)), R(A);
	for (int i = 0; i < n-1; i++) {
		// making vi
		Vector<long double> vi(n-i);
		for (int j = i; j < n; j++)
			vi.at(j-i) = A.elem(j, i);
		vi.at(0) += sign(vi.at(0)) * vi.norm();
		
		// making H_i
		
		// _H_i (R⁽ⁿ⁻ⁱ⁾ˣ⁽ⁿ⁻ⁱ⁾) = I - 2/|vi|² vi viᵀ
		Vector<long double> _H_i(identity(n-i) - (2/vi.norm2()) * vi * vi.transpose());
		
		// H_i (Rⁿˣⁿ)    = (1   0 )
		//                 (0 _H_i)
		Matrix<long double> H_i(n,n);
		for (int r = 0; r < n; r++)
			for (int c = 0; c < n; c++)
				H_i[r][c] = r < i || c < i? (r==c? 1: 0): _H_i[r-i][c-i];
		
		// Updating Q and R
		R = H_i * R;
		Q = Q * H_i;
	}

	return {Q, R};
}