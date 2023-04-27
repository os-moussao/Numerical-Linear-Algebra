#pragma once

#include "../Matrix/Matrix.hpp"

/**
 * Calculating the QR factorization of a matrix A ∈ ℝ𝑛x𝑛, such that: A = product of Q and R
 * where Q is a an orthogonal matrix and R is an upper triangular matrix
 * using the Gram-Schmidt process
 */

// A = [x1, x2, ... xn]
template<typename T>
pair<Matrix<long double>, Matrix<long double>> GramSchmidtQR(Matrix <T> &A) { // QR decomposition using Gram-Schmidt process
	assert(A.rows() && A.rows()==A.cols());

	int n = A.rows();
	Matrix<long double> Q(n,n), R(n,n);

	for (int i = 0; i < n; i++) {
		Vector<long double> xi = A.col(i);
		for (int j = i-1; j >= 0; j--) {
			Vector<long double> qj = Q.col(j);
			Vector<long double> proj_xi_qj = Proj(xi, qj);
			R[j][i] = proj_xi_qj.norm();
			xi = xi - proj_xi_qj;
		}
		R[i][i] = xi.norm();
		xi = xi.unity();
		for (int j = 0; j < n; j++)
			Q[j][i] = xi.at(j);
	}

	return {Q, R};
}