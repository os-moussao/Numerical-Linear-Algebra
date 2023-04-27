#include "GramSchmidtQR.hpp"
#include "HouseholderTransformationsQR.hpp"

int main() {
	int n; cin >> n;
	Matrix<long double> A(n,n);
	cin >> A;
	
	{
		cout << "Gram-Schmidt Process:\n\n";
		
		auto [Q, R] = GramSchmidtQR(A);

		cout << "Q =\n" << Q << endl; // orthogonal matrix
		cout << "R =\n" << R << endl; // upper triangular

		Matrix<long double> A_ = Q*R;
		cout << "QR =\n" << A_ << endl; // should be equal to A
	}

	cout << endl;

	{
		cout << "Househoder Transformations:\n\n";
	
		auto [Q, R] = QR_HouseHolder(A);

		cout << "Q =\n" << Q << endl; // orthogonal matrix
		cout << "R =\n" << R << endl; // upper triangular

		Matrix<long double> A_ = Q*R;
		cout << "QR =\n" << A_ << endl; // should be equal to A

	}
}

/*
example input:
3
1 1 1
1 2 3
1 3 3
*/