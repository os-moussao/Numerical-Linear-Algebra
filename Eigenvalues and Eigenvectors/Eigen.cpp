#include "Eigenvalues.hpp"
#include "PowerMethod.hpp"
#include "ShiftedInversePower.hpp"
#include "Eigen.hpp"

int main() {
	int n; cin >> n;
	Matrix<long double> A(n,n);
	cin >> A;

	auto [eigenvalues, eigenvectors] = Eigen(A);

	for (int i = 0; i < n; i++) {
		long double lam = eigenvalues[i];
		Vector<long double> v = eigenvectors.col(i);

		// assert A v = λ v
		Vector<long double> check(A*v - lam*v);
		for (int j = 0; j < n; j++)
			assert(abs(check.at(j)) <= 1e-10);
		
		cout << "eigenvalue = " << lam << endl;
		cout << "corresponding eigenvector =\n" << v << endl << endl;

	}
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