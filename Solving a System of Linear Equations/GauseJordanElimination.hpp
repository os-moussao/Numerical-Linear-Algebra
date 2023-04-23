#pragma once

#include "../Matrix/Matrix.hpp"

/**
 * Solving a system of linear equations, with $n equations and $n unknowns
 * using gause-jordan elimination algorithm with partial pivoting
 */

template<typename U, typename V>
Matrix<long double> GauseJordanElim(Matrix<U> &factors, Matrix<V> &rhs) {
	assert(factors.rows()==factors.cols() && factors.rows()==rhs.rows() && rhs.cols()==1);
	
	int n = factors.rows(), m = n+1;
	Matrix<long double> sys(n, m);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			sys[i][j] = factors[i][j];
		sys[i][n] = rhs[i][0];
	}

	for (int piv = 0; piv < n; piv++) {
		for (int i = piv+1; i < n; i++) { // partial pivoting, for better precision
			if (abs(sys[piv][piv]) < abs(sys[i][piv])) {
				sys[piv].swap(sys[i]);
			}
		}

		if (sys[piv][piv]==0) {
			// this means either the system has no solution, or infinite solutions
			// and determinant(factors) = 0
			return Matrix<long double>();
		}
		
		if (sys[piv][piv] != 1.) {
			long double sys_piv_piv = sys[piv][piv];
			for (int i = 0; i < m; i++)
				sys[piv][i] /= sys_piv_piv;
		}

		for (int i = 0; i < n; i++) {
			if (i == piv || sys[i][piv] == 0.) continue ;
			long double fact = -sys[i][piv];
			for (int j = 0; j < m; j++)
				sys[i][j] += fact * sys[piv][j];
		}
	}

	Matrix<long double> solution(n,1);
	for (int i = 0; i < n; i++)
		solution[i][0] = sys[i][n];
	return solution;
}