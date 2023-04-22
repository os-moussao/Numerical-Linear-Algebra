#include <bits/stdc++.h>
using namespace std;

template<typename T>
class Matrix {
	int n, m;
	vector<vector<T> > matrix;
	
	public:
	Matrix(int n, int m, T def = 0): n(n), m(m), matrix(n, vector<T>(m, def)){}
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


	T determinant() const { return determinant_naive(); }

	Matrix identity() const {
		assert(n==m);
		Matrix<T> id(n,n);
		for (int i = 0; i < n; i++)
			id[i][i] = 1;
		return id;
	}

	private:
	
	T determinant_naive() const {
		assert(n==m);
		
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
			det += matrix[0][k] * (k&1? -sub.determinant(): sub.determinant());
			r = 0, c = 0;
		}

		return det;
	}

};

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