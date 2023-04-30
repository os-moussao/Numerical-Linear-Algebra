#include "Eigenvalues.hpp"
#include "PowerMethod.hpp"

int main() {
	int n; cin >> n;
	Matrix<long double> A(n,n);
	cin >> A;

	vector<long double> eigvalues = eigenvalues(A);
	sort(eigvalues.begin(), eigvalues.end(), [&] (auto &x, auto y) {
		return abs(x) > abs(y);
	});

	for (long double x: eigvalues)
		cout << x << '\n';

	cout << "\nPower Method:" << endl;

	auto[lam, firstEigenvector] = PowerMethod(A);

	cout << "\neigenvalue with largest magnitude = " << lam << endl; // should be equal to eigvalues[0] 
	cout << "\nx_0 =\n" << firstEigenvector << endl;

	Vector<long double> a(A*firstEigenvector), b(eigvalues[0] * firstEigenvector);

	// a,b should be equal
	cout << "A x_0 =\n" << a << endl;
	cout << "lam x_0 =\n" << b << endl;
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