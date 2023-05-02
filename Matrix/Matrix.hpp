#pragma once

#include <bits/stdc++.h>
using namespace std;

template<typename T> class Matrix;
template<typename T> ostream &operator<< (ostream &os, Matrix<T> &m);

#define Vector Matrix

template<typename T>
class Matrix {
	
	protected:
		int n, m;
		vector<vector<T> > matrix;
	
	public:
	Matrix(int n = 0, int m = 1, T def = 0): n(n), m(m), matrix(n, vector<T>(m, def)){}
	template<typename V> Matrix(const Matrix<V> &cpy) { *this = cpy; }

	int rows() const { return n; }
	int cols() const { return m; }
	vector<T> &operator [] (int i) { return matrix[i]; }
	T elem(unsigned int i, unsigned int j = 0) const {
		assert(i < (const unsigned int)n && j < (const unsigned int)m);
		return matrix[i][j];
	}

	Matrix &operator = (const Matrix &cpy) {
		this->n = cpy.n;
		this->m = cpy.m;
		this->matrix = cpy.matrix;
		return *this;
	}

	template<typename X>
	Matrix &operator=(const Matrix<X> &b) {
		if (n != b.rows() || m != b.cols())
			n=b.rows(), m=b.cols(), matrix.resize(b.rows(), vector<T>(b.cols()));
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				matrix[i][j] = b.elem(i,j);
		return *this;
	}

	bool operator==(const Matrix<T> &B) const {
		if (n != B.rows() || m != B.cols())
			return false;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				if (matrix[i][j] != B.elem(i, j))
					return false;
			}
		}
		return true;
	}

	Matrix operator *(const Matrix<T> &right) const {
		assert(cols() == right.rows());
		Matrix<T> mult(rows(), right.cols());
		for (int i = 0; i < mult.rows(); i++)
			for (int j = 0; j < mult.cols(); j++)
				for (int k = 0; k < cols(); k++)
					mult[i][j] += matrix[i][k] * right.elem(k, j);
		return mult;
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

	Matrix operator + (const Matrix<T> &B) const {
		assert(n==B.rows() && m==B.cols());
		Matrix<T> A(n, m);
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				A[i][j] = matrix[i][j] + B.elem(i, j);
		return A;
	}

	Matrix operator - (const Matrix<T> &B) const {
		assert(n==B.rows() && m==B.cols());
		Matrix<T> A(n, m);
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				A[i][j] = matrix[i][j] - B.elem(i,j);
		return A;
	}

	template<typename V>
	friend Matrix<V> operator *(long double scalar, const Matrix<V> &A) {
		Matrix<V> B(A);
		for (int i = 0; i < B.rows(); i++)
			for (int j = 0; j < B.cols(); j++)
				B[i][j] *= scalar;
		return B;
	}

	T determinant() const {
		assert(n==m);
		return determinant_gauss_elimination();
	}

	Matrix identity() const {
		assert(n==m);
		Matrix<T> id(n,n);
		for (int i = 0; i < n; i++)
			id[i][i] = 1;
		return id;
	}

	Matrix transpose() const {
		Matrix<T> mat(m, n);
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				mat[j][i] = matrix[i][j];
		return mat;
	}

	Matrix inverse() {
		assert(n==m);

		Matrix<T> inverse = identity();
		Matrix<T> mat(*this);

		for (int piv = 0; piv < n; piv++) {
			int optimal_pivot = piv;
			for (int i = piv+1; i < n; i++) { // partial pivoting, for better precision
				if (abs(mat[optimal_pivot][piv]) < abs(mat[i][piv])) {
					optimal_pivot = i;
				}
			}
			mat[piv].swap(mat[optimal_pivot]);
			inverse[piv].swap(inverse[optimal_pivot]);


			if (abs(mat[piv][piv]) < 1e-10) {
				// then the matrix is singular
				return Matrix<T>();
			}
			
			if (mat[piv][piv] != 1) {
				long double mat_piv_piv = mat[piv][piv];
				for (int i = 0; i < m; i++)
					mat[piv][i] /= mat_piv_piv, inverse[piv][i] /= mat_piv_piv;
			}

			for (int i = 0; i < n; i++) {
				if (i == piv || mat[i][piv] == 0.) continue ;
				long double fact = -mat[i][piv];
				for (int j = 0; j < m; j++)
					mat[i][j] += fact * mat[piv][j], inverse[i][j] += fact * inverse[piv][j];
			}
		}

		return inverse;
	}

	// Vector utils

	private:
		void __vec__() const { assert(n && m==1); }

	public:

	int size() const { return __vec__(), n; }

	long double norm2() const { // norm^2 to avoid precision lost by sqrtl
		return __vec__(), dot(*this, *this);
	}
	
	long double norm() const {
		return __vec__(), sqrtl(norm2());
	}

	Vector<long double> unity() const {
		__vec__();
		Vector<long double> u; u = *this;
		long double nrm = norm();
		for (int i = 0; i < n; i++)
			u.at(i) /= nrm;
		return u;
	}

	T &at(unsigned int i) {
		return __vec__(), assert(i < n), this->matrix[i][0];
	}

	Vector<T> col(unsigned int i) const {
		assert(i < m);
		Vector<T> col_i(n);
		for (int r = 0; r < n; r++)
			col_i.at(r) = matrix[r][i];
		return col_i;
	}

	template<typename U, typename V>
	friend long double dot(const Vector<U> &u, const Vector<V> &v) {
		assert(u.size()==v.size());
		long double dotProd = 0;
		for (int i = 0; i < v.size(); i++) {
			dotProd += v.elem(i) * u.elem(i);
		}
		return dotProd;
	}

	template<typename U, typename V>
	friend Vector<long double> Proj(const Vector<U> &a, const Vector<V> &v) {
		assert(a.size()==v.size());
		Vector<long double> u = v.unity();
		return dot(a,u) * u;
	}

	private:

	T determinant_gauss_elimination() const {
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

// gauss elimination for integer type
template<>
int Matrix<int>::determinant_gauss_elimination() const {
	Matrix<int> mat(*this);
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

Matrix<long double> identity(int n) {
	Matrix<long double> id(n,n);
	for (int i = 0; i < n; i++)
		id[i][i] = 1;
	return id;
}


// IO utils
template<typename T>
ostream &operator<< (ostream &os, const Matrix<T> &m) {
	os << '[';
	for (int i = 0; i < m.rows(); i++) {
		if (i) os << " ";
		os << '[';
		for (int j = 0; j < m.cols(); j++) {
			if (j) os << ", ";
			os << m.elem(i, j);
		}
		os << ']';
		if (i < m.rows()-1) os << ",\n";
	}
	return os << ']';
}

template<typename T>
ostream &operator<< (ostream &os, Matrix<T> &m) {
	return os << (const Matrix<T> &)m;
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