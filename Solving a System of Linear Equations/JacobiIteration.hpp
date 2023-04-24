#pragma once

#include "../Matrix/Matrix.hpp"

/**
 * Solving a system of linear equations using Jacobi Iteration
 * This is meant to solve for large system of equations, where O(N^3) (using direct solvers) is too slow
 * Better to use iterative methods on diagonally dominant matrices, as they're guaranteed to converge
 */

#define EPS 1e-10

template<typename T>
bool DiagonallyDominant(Matrix<T> &matrix) {
	assert(matrix.rows()==matrix.cols());
	bool ok = 0;
	for (int i = 0; i < matrix.rows(); i++) {
		long double diff = abs(matrix[i][i]);
		for (int j = 0; j < matrix.cols(); j++) if (j!=i)
			diff -= abs(matrix[i][j]);
		if (diff < 0.) return false;
		if (diff > 0.) ok=1;
	}
	return ok;
}

template<typename U, typename V>
long double _error(Matrix<U> &factors, Matrix<V> &rhs, Matrix<long double> &solution) {
	long double max_error = -10;
	for (int i = 0; i < factors.rows(); i++) {
		long double res = 0.;
		for (int j = 0; j < factors.cols(); j++)
			res += factors[i][j] * solution[j][0];
		max_error = max(max_error, abs((long double)rhs[i][0] - res));
	}
	return max_error;
}

template<typename U, typename V>
Matrix<long double> Jacobi(Matrix<U> &factors, Matrix<V> &rhs, long double tolerance = EPS) {
	assert(factors.rows()==factors.cols() && factors.rows()==rhs.rows() && rhs.cols()==1);

	int n = factors.rows();
	int max_iterations = 1e7/(2*n*n + n);
	// complexity of each iteration is O(N^2)
	// we solve for [ max_iterations * O(N^2) = 10^7 instructions (~ 1 second in modern computers) ]

	Matrix<long double> solution(n,1,0.); // initial guess of 0's

	for (int k = 0; k < max_iterations && _error(factors, rhs, solution) > tolerance; k++) {
		Matrix<long double> new_solution(n,1);
		for (int i = 0; i < n; i++) {
			new_solution[i][0] = rhs[i][0];
			for (int j = 0; j < n; j++) if (j != i)
				new_solution[i][0] -= factors[i][j] * solution[j][0];
			new_solution[i][0] /= factors[i][i];
		}
		solution = new_solution;
	}

	return solution;
}