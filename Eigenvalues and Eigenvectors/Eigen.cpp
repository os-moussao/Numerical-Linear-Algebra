#include "Eigenvalues.hpp"
#include "PowerMethod.hpp"
#include "ShiftedInversePower.hpp"

int main() {
	int n; cin >> n;
	Matrix<long double> A(n,n);
	cin >> A;

	vector<long double> eigvalues = eigenvalues(A);

	// eigenvalue, eigenvector pairs
	vector<pair<long double, Vector<long double>>> eigen(n);
	for (int i = 0; i < n; i++) {
		// get the closest eigenvalue to eigvalues[i] and the corresponding eigenvector
		eigen[i] = ShiftedInversePower(A, eigvalues[i]-0.01);
	}

	for (auto &[lambda, v]: eigen) {
		// assert A v = λ v
		Vector<long double> Av(A * v);
		for (int i = 0; i < n; i++)
			assert(abs(Av.at(i) - lambda * v.at(i)) <= 1e-10);

		cout << "eigenvalue = " << lambda << endl;
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