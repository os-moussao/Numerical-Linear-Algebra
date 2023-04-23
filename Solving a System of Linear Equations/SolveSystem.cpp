#include "GauseJordanElimination.hpp"

int main() {
	int n; cin >> n;
	Matrix<long double> a(n,n), b(n, 1);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			cin >> a[i][j];
		cin >> b[i][0];
	}

	Matrix<long double> solution = GauseJordanElim(a, b);
	if (!solution.rows()) {
		cout << "The system either has 0 or infinite solutions\n";
		return 1;
	}
	cout << solution << endl;
}

/*
	input example:
3
1      1  1 2
0      1  1 1
0.0001 1  2 1.9999
*/