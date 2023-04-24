#include "Matrix.hpp"

int main() {
	int n; cin >> n;
	Matrix<int> a(n, n);
	cin >> a;

	cout << "a =\n" << a << endl << endl;
	cout << "determinant of a = " << a.determinant() << endl << endl;

	Matrix<long double> b; b = a;
	cout << "b =\n" << b << endl << endl;
	cout << "determinant of b = " << b.determinant() << endl;

	Matrix<long double> ib = b.inverse();
	cout << "\nInverse of b =\n";
	if (ib.rows()) {
		Matrix<long double> id = b.identity();
		Matrix<long double> id1 = ib * b, id2 = b * ib;
		
		// assert solution is very close
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				assert((abs(id[i][j]-id1[i][j]) <= 1e-6) && (abs(id[i][j]-id2[i][j]) <= 1e-6));
			}
		}
		cout << ib << endl;
	} else {
		cout << "b is a singular matrix !\n";
	}
}