#pragma once

#include <bits/stdc++.h>
using namespace std;

template<typename T> class Matrix;
template<typename T> ostream &operator<< (ostream &os, Matrix<T> &m);

template<typename T>
class Matrix {
	int n, m;
	vector<vector<T> > matrix;
	
	public:
	Matrix(int n = 0, int m = 0, T def = 0): n(n), m(m), matrix(n, vector<T>(m, def)){}
	Matrix(const Matrix<T> &cpy) { *this = cpy; }

	int rows() const { return n; }
	int cols() const { return m; }
	vector<T> &operator [] (int i) { return matrix[i]; }

	Matrix &operator = (const Matrix &cpy) {
		this->n = cpy.n;
		this->m = cpy.m;
		this->matrix = cpy.matrix;
		return *this;
	}

	template<typename X>
	Matrix &operator=(Matrix<X> &b) {
		if (n != b.rows() || m != b.cols())
			n=b.rows(), m=b.cols(), matrix.resize(b.rows(), vector<T>(b.cols()));
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				matrix[i][j] = b[i][j];
		return *this;
	}

	Matrix operator *(const Matrix<T> &right) const {
		assert(cols() == right.rows());
		Matrix<T> mult(rows(), right.cols());
		for (int i = 0; i < mult.rows(); i++)
			for (int j = 0; j < mult.cols(); j++)
				for (int k = 0; k < cols(); k++)
					mult[i][j] += matrix[i][k] * right[k][j];
	}

	Matrix operator ^ (unsigned int exp) const {
		assert(n==m);
		Matrix<T> pr = identity();
		Matrix<T> mat = *this;
		while (exp) {
			if (exp&1) pr = pr * mat;
			mat = mat * mat;
			exp >>= 1;
		}
		return pr;
	}


	T determinant() {
		assert(n==m);
		return determinant_gause_elimination();
	}

	Matrix identity() const {
		assert(n==m);
		Matrix<T> id(n,n);
		for (int i = 0; i < n; i++)
			id[i][i] = 1;
		return id;
	}

	private:

	T determinant_gause_elimination() {
		// converting to long double
		Matrix<long double> mat;
		mat = *this;

		long double det = 1;
		int sign = 1;

		for (int piv = 0; piv < n; piv++) {
			for (int i = piv+1; i < n; i++) { // partial pivoting
				if (abs(mat[piv][piv]) < abs(mat[i][piv])) {
					mat[piv].swap(mat[i]);
					sign *= -1;
					break ;
				}
			}
			if (mat[piv][piv]==0) // then determinant is 0
				return 0;
			for (int i = piv+1; i < n; i++) {
				if (mat[i][piv] == 0.) continue;
				long double fact = mat[i][piv] / mat[piv][piv];
				for (int j = 0; j < n; j++) {
					mat[i][j] = mat[i][j] - fact * mat[piv][j];
				}
			}
		}
		for (int i = 0; i < n; i++)
			det *= mat[i][i];
		return det * sign;
	}
	
	T determinant_naive() const {
		if (n==1)
			return matrix[0][0];
		
		T det = 0;
		int r = 0, c = 0;
		Matrix<T> sub(n-1,n-1);
		for (int k = 0; k < n; k++) {
			for (int i = 1; i < n; i++) {
				for (int j = 0; j < n; j++) {
					if (j == k) continue;
					sub[r][c] = matrix[i][j];
					if (++c == n-1)
						r++, c = 0;
				}
			}
			det += matrix[0][k] * (k&1? -sub.determinant_naive(): sub.determinant_naive());
			r = 0, c = 0;
		}

		return det;
	}

};

// gause elimination for integer type
template<>
int Matrix<int>::determinant_gause_elimination() {
	Matrix<long long> mat; mat = *this;
	long long det = 1;
	long long fact = 1;
	for (int piv = 0; piv < n; piv++) {
		if (!mat[piv][piv]) {
			for (int i = piv+1; i < n; i++) {
				if (mat[i][piv]) {
					mat[piv].swap(mat[i]);
					fact *= -1;
					break ;
				}
			}
			if (!mat[piv][piv])
				return 0;
		}
		for (int i = piv+1; i < n; i++) {
			if (!mat[i][piv]) continue;
			int lc = lcm(abs(mat[i][piv]), abs(mat[piv][piv]));
			int fact_i = lc / mat[i][piv], fact_piv = lc / mat[piv][piv];
			for (int j = 0; j < n; j++)
				mat[i][j] = mat[i][j] * fact_i - mat[piv][j] * fact_piv;
			fact *= fact_i;
		}
	}
	for (int i = 0; i < n; i++)
		det *= mat[i][i];
	return det /= fact, det;
}

template<typename T>
ostream &operator<< (ostream &os, Matrix<T> &m) {
	os << '[';
	for (int i = 0; i < m.rows(); i++) {
		if (i) os << " ";
		os << '[';
		for (int j = 0; j < m.cols(); j++) {
			if (j) os << ", ";
			os << m[i][j];
		}
		os << ']';
		if (i < m.rows()-1) os << ",\n";
	}
	return os << ']';
}

template<typename T>
istream &operator>> (istream &is, Matrix<T> &m) {
	for (int i = 0; i < m.rows(); i++) {
		for (int j = 0; j < m.cols(); j++) {
			is >> m[i][j];
		}
	}
	return is;
}