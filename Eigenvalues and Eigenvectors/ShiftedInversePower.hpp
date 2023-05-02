#pragma once

#include "PowerMethod.hpp"

/**
 * Shifted Inverse Power Method for finding the eigenvalue closest to some α, and its coresponding eigenvector.
 * 
 * if λ is an eigenvalue of A, with A v = λ v
 * then 1/λ is an eigenvalue of inverse(A) and inverse(A) = 1/λ v
 * 
 * so if applying the Power Method on inverse(A) gives the eigenvalue μ of inverse(A) with largest magnitude
 * then λ = 1/μ is the eigenvalue with smallest magnitude of A
 * 
 * also if A v = λ v, then for any α: (A - α I) v = (λ - α) v
 * 
 * let α be some constant:
 * applying Power Method on inverse(A - α I) give us the eigenvalue ω
 * then 1/ω = λ - α is the smallest eigenvalue of (A - α I) in magnitude
 * therefore λ = 1/ω + α is the closest eigenvalue of A to α.
 */

template<typename T>
pair<long double, Vector<long double>> ShiftedInversePower(const Matrix<T> &A, long double alpha) {
	assert(A.rows() && A.rows()==A.cols());

	// creating B = inverse(A - alpha I)
	Matrix<long double> B; B = A;
	for (int i = 0; i < A.rows(); i++)
		B[i][i] -= alpha;
	B = B.inverse();

	auto [omega, eigenVector] = PowerMethod(B);

	return {1/omega + alpha, eigenVector};
}