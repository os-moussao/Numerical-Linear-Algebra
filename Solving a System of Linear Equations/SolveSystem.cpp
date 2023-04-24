#include "GaussJordanElimination.hpp"
#include "JacobiIteration.hpp"
#include "GaussSeidelIteration.hpp"

int main() {
	int n; cin >> n;
	Matrix<long double> a(n,n), b(n, 1);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			cin >> a[i][j];
		cin >> b[i][0];
	}

	Matrix<long double> gaussSolution = GaussJordanElim(a, b);
	if (!gaussSolution.rows())
		cout << "Gauss Jordan Elimination: The system either has 0 or infinite solutions\n";
	else
		cout << "Gauss Jordan Elinination:\n" << gaussSolution << endl;

	Matrix<long double> jacobiSolution = Jacobi(a, b);
	cout << "Jordan Iteration:\n" << jacobiSolution << endl;

	Matrix<long double> gaussSeidleSolution = GaussSeidel(a, b);
	cout << "Gauss-Seidle Iteration:\n" << jacobiSolution << endl;

}

/*
	input example:
3
1      1  1 2
0      1  1 1
0.0001 1  2 1.9999
*/