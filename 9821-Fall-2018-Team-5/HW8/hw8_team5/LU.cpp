#include "header.hpp"

/*(2) Write C++ code called lu no pivoting and lu row pivoting to compute the LU decom
position without pivoting of a matrix and the LU decomposition with row pivoting of
a matrix. The row pivoted LU should satisfy P A = LU.
*/

std::tuple <mat, mat > lu_no_pivoting(mat A)
{
	mat L = mat::Identity(A.rows(), A.cols());
	mat U = mat::Zero(A.rows(), A.cols());
	auto n = A.rows();

	for (auto i = 0; i < n - 1; i++)
	{
		for (auto k = i; k < n; k++)
		{
			U(i, k) = A(i, k);
			L(k, i) = A(k, i) / U(i, i);
		}
		for (auto j = i + 1; j < n; j++)
		{
			for (auto k = i + 1; k < n; k++)
			{
				A(j, k) = A(j, k) - L(j, i)*U(i, k);
			}
		}
	}
	
	U(n - 1, n - 1) = A(n - 1, n - 1);
	return std::make_tuple(L, U);
}


std::tuple <permutation, mat, mat> lu_row_pivoting(mat A)
{
	Eigen::PartialPivLU<mat> lu(A);
	//generate the L matrix
	mat L = mat::Identity(A.rows(), A.cols());
	L.triangularView < Eigen::StrictlyLower > () = lu.matrixLU();
	//generate the U matrix
	mat U = lu.matrixLU().triangularView<Eigen::Upper>();
	//generate the permutation matrix
	permutation P = lu.permutationP();
	return std::make_tuple(P, L, U);
}