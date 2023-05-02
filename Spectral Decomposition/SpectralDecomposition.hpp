#pragma once 

#include "../Eigenvalues and Eigenvectors/Eigen.hpp"

/**
 * Spectral Decomposition for Symmetric Matrices
 * 
 * Important properties of symmetric matrices:
 * 1. They have real eigenvalues
 * 2. Their eigenvectors corresponding to different eigenvalues are orthogonal
 * 
 * Let A (in Rⁿˣⁿ) be symmetric:
 * Let Q = [q_1 | q_2 | ... | q_n], be the eigenvecors of A
 * Let V = diag(λ_1, λ_2, ... , λ_n) be the eigenvalues of A in a diagonal
 * A can be formed ad A = Q V Qᵀ
 */

struct Spect {
	Matrix<long double> Q, V, QT;
};

template<typename T>
Spect SpectralDecomp(const Matrix<T> &A) {
	// asserting symmetricity 
	assert(A.rows() && (A == A.transpose()));
	
	int n = A.rows();

	auto [evalues, Q] = Eigen(A);
	
	Matrix<long double> V(n, n);
	for (int i = 0; i < n; i++)
		V[i][i] = evalues[i];
	
	return (Spect){Q, V, Q.transpose()};
}
