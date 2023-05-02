#include "SpectralDecomposition.hpp"

int main() {
	int n; cin >> n;
	Matrix<long double> A(n, n);
	cin >> A;

	// A has to be symmetric
	assert(n && (A == A.transpose()));

	Spect decomp = SpectralDecomp(A);

	Matrix<long double> &Q = decomp.Q;
	Matrix<long double> &V = decomp.V;
	Matrix<long double> &QT = decomp.QT;

	Matrix<long double> check(Q*V*QT - A);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			assert(abs(check[i][j]) <= 1e-10);

	cout << "Q =\n" << Q << "\n\n";
	cout << "V =\n" << V << "\n\n";
	cout << "Qᵀ =\n" << QT << "\n\n";

	cout << "Q * V * Qᵀ =\n" << Q*V*QT << endl;
}