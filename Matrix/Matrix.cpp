#include "Matrix.hpp"

int main() {
	int n; cin >> n;
	Matrix<int> a(n, n);
	cin >> a;

	cout << "a =\n" << a << endl << endl;
	cout << "determinant of a = " << a.determinant() << endl;

	Matrix<double> b(n, n);
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			b[i][j] = a[i][j];
	cout << "determinant of b = " << b.determinant() << endl;
}