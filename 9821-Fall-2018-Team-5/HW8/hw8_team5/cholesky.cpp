//cholesky
/*(3) Write C++ code for the Cholesky decomposition of an spd matrix. Return U such
that UtU = A.
*/

#include "header.hpp"
mat cholesky(mat A)
{
	int n = A.rows();
	mat U = mat::Zero(A.rows(), A.cols());
	for (auto i = 0; i < (n - 1); i++)
	{
		U(i, i) = sqrt(A(i, i));
		for (auto k = i + 1; k < n; k++)
		{
			U(i, k) = A(i, k) / U(i, i);
		}
		for (auto j = i + 1; j < n; j++)
			for (auto m = j; m < n; m++)
				A(j, m) = A(j, m) - U(i, j)*U(i, m);
	}
	U(n-1, n-1) = sqrt(A(n-1, n-1));
	return U;
}

vec cholesky_linear(mat &A, vec &b)
{
	mat U(A.rows(), A.cols());
	U = cholesky(A);
	vec y(A.rows());
	y = forward_subst(U.transpose(), b);
	vec x(A.rows());
	x = backward_subst(U, y);
	return x;
}
