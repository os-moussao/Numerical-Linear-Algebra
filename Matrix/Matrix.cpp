#include "Matrix.hpp"

int main() {
	int n; cin >> n;
	Matrix<int> a(n, n);
	cin >> a;

	cout << "a =\n" << a << endl << endl;
	cout << "determinant of a = " << a.determinant() << endl << endl;

	Matrix<double> b; b = a;
	cout << "b =\n" << b << endl << endl;
	cout << "determinant of b = " << b.determinant() << endl;
}