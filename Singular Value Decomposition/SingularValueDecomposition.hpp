#pragma once

#include "../Eigenvalues and Eigenvectors/Eigen.hpp"

struct Svd {
	Matrix<long double> U, Sigma, VT; // U Σ Vᵀ
};

template<typename T>
Svd SVD(const Matrix<T> &A) {
	assert(A.rows() && A.cols());

	int n = A.rows(), m = A.cols();

	// Getting eigenvalues and eigenvectors of AᵀA
	auto [svalues, V] = Eigen(A.transpose() * A);

	// getting singular values of A (σᵢ = √λᵢ)
	for (long double &lambda: svalues)
		lambda = sqrtl(lambda);

	// Computing Σ
	Matrix<long double> Sigma;
	Sigma = Matrix<long double>(n, m);
	for (int i = 0; i < min(n,m); i++) {
		Sigma[i][i] = svalues[i];
	}

	for (int i = min(n, m); i < m; i++) {
		for (int j = 0; j < m; j++) {
			V[j][i] = 0;
		}
	}

	// Computing U, u_i = 1/σ_i A v_i
	Matrix<long double> U(n, n);
	for (int i = 0; i < n; i++) {
		Vector<long double> u_i = i < min(n, m)? // in range of diagonal
									(1/Sigma[i][i]) * A * V.col(i):
									Vector<long double>(n);
		for (int j = 0; j < n; j++)
			U[j][i] = u_i.at(j);
	}

	return (Svd) {U, Sigma, V.transpose()};
}


void isOrth(const Matrix<long double> &V, string name) {
	assert(V.cols()==V.rows());
	int n = V.cols();
	cerr << "is " << name << " orthogonal? ";
	for (int i = 0; i < n; i++) {
		Vector<long double> e_i = V.col(i);
		for (int j = i+1; j < n; j++) {
			Vector<long double> e_j = V.col(j);
			assert(abs(dot(e_i,e_j)) <= 1e-10);
		}
	}
	cerr << "Yes\n";
}