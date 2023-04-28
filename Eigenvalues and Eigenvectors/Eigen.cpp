#include "Eigenvalues.hpp"

int main() {
	int n; cin >> n;
	Matrix<long double> A(n,n);
	cin >> A;

	vector<long double> eigen = eigenvalues(A);
	for (long double x: eigen)
		cout << x << '\n';
}
/*
example input:

symetrix matrix:
3
1 1 1
1 2 3
1 3 3

asymetric matrix:
3
1 2 1
1 2 3
1 3 3

online eigenvalues calcuator: https://matrixcalc.org/vectors.html
*/