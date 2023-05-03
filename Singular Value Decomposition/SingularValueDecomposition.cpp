#include "SingularValueDecomposition.hpp"

int main() {
	int n, m; cin >> n >> m;
	Matrix<long double> A(n,m);
	cin >> A;

	Svd decomp = SVD(A);

	Matrix<long double> &U = decomp.U, &Sigma = decomp.Sigma, &VT = decomp.VT;


	cout << "U =\n" << U << endl;
	cout << "Σ =\n" << Sigma << endl;
	cout << "Vᵀ =\n" << VT << endl;

	cout << "U Σ Vᵀ =\n";
	cout << U * Sigma * VT << endl; // Should be equal to A

	isOrth(U, "U"); // 
	isOrth(VT.transpose(), "V");
}

/*
example input:

3 5
1 2 3 14 -3
7 10 1 23 29
-1 5 3 8 26


5 3
1 2 3
7 10 1
-1 5 3
14 -3 8
23 29 1

*/