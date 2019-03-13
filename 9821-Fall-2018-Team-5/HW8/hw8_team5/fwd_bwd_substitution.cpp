#include "header.hpp"

/*  1
Write C++ codes for backward and forward substitutions, called forward subst and
backward subst.*/

vec forward_subst(const mat & L, const vec & b)
{
	auto x = b;
	x(0) = b(0) / L(0, 0);

	for (auto j = 1; j < b.size(); j++)
	{
		double sum = 0;
		for (auto k = 0; k < j; k++)
		{
			sum = sum + L(j, k) * x(k);
		}
		x(j) = (b(j) - sum) / L(j,j);
	}
	return x;
}

vec backward_subst(const mat & U, const vec & b)
{
	auto x = b;
	auto n = b.size();
	x(n - 1) = b(n - 1) / U(n - 1, n - 1);

	for (auto j = n - 2; j >= 0; j--)
	{
		double sum = 0;
		for (auto k = (j + 1); k < n; k++)
		{
			sum = sum + U(j, k) * x(k);
		}
		x(j) = (b(j) - sum) / U(j, j);
	}
	return x;
}