#include "Matrix.hpp"

int main() {
	int n, m; cin >> n >> m;
	Matrix<double> a(n, n);
	cin >> a;
	cout << "a =\n" << a << endl << endl;
	cout << "determinant of a = " << a.determinant() << endl;
}