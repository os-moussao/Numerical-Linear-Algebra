#pragma once

#include "../Matrix/Matrix.hpp"

/**
 * Power method for finding the eigenvalue with the largest magnitude, and its corresponding eigenvector
 * 
 * if A ∈ ℝ𝑛x𝑛 has n independent eigenvectors, A x_i = λ_i x_i (for convenience, let |λ_i| > |λ_i+1|)
 * then any vector v ∈ ℝ𝑛 can be written as a linear combination of all eigenvectors of A
 * v = c_0 x_0 + c_1 x_1 + ... + c_n-1 x_n-1
 * v_k = A^k v, converges to λ_0^k c_0 x_0, A^k v = λ_0^k c_0 x_0
 * so v_k is in the direction of x_0
 * 
 * first we only care about direction of x_0, so to avoid numerical errors
 * after each iteration we will normalize v_k
 * but before doing that we keep track of λ_0 which is equal to dotProd(v_k+1, v_k)
 * v_k+1 = A v_k
 * λ_0 = dot(v_k+1, v_k)
 * v_k+1 /= |v_k+1|
 * 
 * we stop when we exceed maximum iterations, or when we completely converge (λ_0 / v_k do not change)
 */

template<typename T>
pair<long double, Vector<long double>> PowerMethod(const Matrix<T> &A) {
	assert(A.rows() && A.rows()==A.cols());

	int n = A.rows();

	int max_iterations = 10000;

	long double lambda = 0;
	Vector<long double> v(n,1,1);
	v = v.unity();

	for (int k = 0; k < max_iterations; k++) {
		Vector<long double> w = A * v;

		long double new_lambda = dot(w, v);
		
		if (new_lambda == lambda) // convergence
			break ;
		
		v = w.unity(), lambda = new_lambda;
	}

	return {lambda, v};
}