#pragma once

#include "Eigenvalues.hpp"
#include "ShiftedInversePower.hpp"

/**
 * Get eigenvalues and a Matrix of eigenvectors of A 
 * Sorted in decreasing order of magnitudes
 */

template<typename T>
pair<vector<long double>, Matrix<long double>> Eigen(const Matrix<T> &A) {
	assert(A.rows() && A.rows()==A.cols());

	int n = A.rows();
	vector<long double> eigvalues = eigenvalues(A);

	sort(eigvalues.begin(), eigvalues.end(), [&](auto &x, auto &y) {
		return abs(x) > abs(y);
	});

	Matrix<long double> eigvectors(n, n);
	for (int i = 0; i < n; i++) {
		Vector<long double> v = ShiftedInversePower(A, eigvalues[i]-0.0001).second;
		for (int j = 0; j < n; j++)
			eigvectors[j][i] = v.at(j);
	}

	return {eigvalues, eigvectors};
}
