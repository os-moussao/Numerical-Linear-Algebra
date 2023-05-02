#pragma once

#include "../Matrix/Matrix.hpp"
#include "../QR factorization/HouseholderTransformationsQR.hpp"

#ifndef EPS
#define EPS 1e-10
#endif

/**
 * Finding the approximate eigen values of A ∈ ℝ𝑛x𝑛 using QR decomposition
 * A_0 = A = Q_0 R_0
 * for k = 0 ... untill we exceed max iterations or A_k is approximate to an upper triangular matrix:
 * 	    A_k+1 = R_k Q_k // A_k+1 is similar to A_k, since A_k+1 = inv(Q_k) A_k Q_k (the similarity transform)
 * the approximate eigen values of A are the diagonal of A_k
 */

template<typename T>
bool isUpperTriangular(Matrix<T> &A, int n, long double tolerance) {
	for (int i = 1; i < n; i++) {
		for (int j = 0; j < i; j++)
			if (abs(A[i][j]) > EPS)
				return false;
	}
	return true;
}

template<typename T>
vector<long double> eigenvalues(const Matrix<T> &A, long double tolerance=EPS) {
	assert(A.rows() && A.rows()==A.cols());

	int n = A.rows();

	int max_iterations = 10000;

	Matrix<long double> A_k = A;

	for (int k = 1; k < max_iterations && !isUpperTriangular(A_k,n,tolerance); k++) {
		auto [Q, R] = QR_HouseHolder(A_k); // A_k = Q_k R_k
		A_k = R * Q; // A_k+1 = R_k * Q_k // A_k <- A_k+1
	}

	vector<long double> eigen(n);
	for (int i = 0; i < n; i++)
		eigen[i] = A_k[i][i];
	
	return eigen;
}