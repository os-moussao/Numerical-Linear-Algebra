#pragma once

#include "JacobiIteration.hpp" // for _error calculation

/**
 * Gauss-Seidel Iteration is an improvement of Jacobi Iteration method
 * where: x_i(k) = (b_i - Σ j:(1->n, j != i) x_j(k-1) * a_i_j) / a_i_i
 * becomes: x_i(k) = (b_i - Σ j:(1->i-1) x_j(k) * a_i_j - Σ j:(i+1->n) x_j(k-1) * a_i_j) / a_i_i
 */

template<typename U, typename V>
Matrix<long double> GaussSeidel(Matrix<U> &factors, Matrix<V> &rhs, long double tolerance = EPS) {
	assert(factors.rows()==factors.cols() && factors.rows()==rhs.rows() && rhs.cols()==1);

	int n = factors.rows();
	int max_iterations = 1e7/(2*n*n + n);
	// complexity of each iteration is O(N^2)
	// we solve for [ max_iterations * O(N^2) = 10^7 instructions (~ 1 second in modern computers) ]

	Matrix<long double> solution(n,1,0.); // initial guess of 0's

	for (int k = 0; k < max_iterations && _error(factors, rhs, solution) > tolerance; k++) {
		Matrix<long double> new_solution(solution);
		for (int i = 0; i < n; i++) {
			new_solution[i][0] = rhs[i][0];
			for (int j = 0; j < n; j++) if (j != i)
				new_solution[i][0] -= factors[i][j] * new_solution[j][0];
			new_solution[i][0] /= factors[i][i];
		}
		solution = new_solution;
	}

	return solution;
}