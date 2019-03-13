#include "header.hpp"

std::tuple<mat, mat, mat> decompose_LDU(const mat&A)
{
	mat L = mat::Zero(A.rows(),A.cols());
	mat D = mat::Zero(A.rows(), A.cols());
	mat U = mat::Zero(A.rows(), A.cols());
	for (auto i = 0; i < A.cols(); i++)
		for (auto j = 0; j < A.cols(); j++) {
			if (i == j)
				D(i, i) = A(i, i);
			else if (i < j)
				U(i, j) = A(i, j);
			else
				L(i, j) = A(i, j);
		}
	return std::make_tuple(L, D, U);
}


std::tuple <vec, int,double>
jacobi(const mat & A, const vec & b, const vec & x_0,
	const double tolerance,
	const StoppingCriterion criterion)
{
	vec x = x_0;
	auto ic = 0;
	auto LDU = decompose_LDU(A);
	auto L = std::get<0>(LDU);
	auto D = std::get<1>(LDU);
	auto U = std::get<2>(LDU);
	auto b_new = D.inverse()*b;
	double res = 0;
	if (criterion == consecutive)
	{
		vec x_old = x_0;
		vec diff = x_0;
		diff(0) = 10;
		while ((diff).norm() >= tolerance)
		{
			x = -D.inverse() * (L*x_old + U*x_old) + b_new;
			diff = x - x_old;
			x_old = x;
			ic = ic + 1;
		}
		res = diff.norm();
		std::cout << std::endl << "norm" << diff.norm() << std::endl;
	}

	else if (criterion == residual)
	{
		vec r_0 = b - A*x_0;
		vec r = r_0;
		while (r.norm() >= tolerance)
		{
			x = -D.inverse() * (L*x + U*x) + b_new;
			r = b - A * x;
			ic = ic + 1;
		}
	}

	return std::make_tuple(x, ic, res);
}



std::tuple <vec, int,double >
jacobi(const mat & A, const vec & b,
	const double tolerance,
	const StoppingCriterion criterion)
{
	vec x_0 = vec::Zero(b.size());
	vec x = x_0;
	auto ic = 0;
	auto LDU = decompose_LDU(A);
	auto L = std::get<0>(LDU);
	auto D = std::get<1>(LDU);
	auto U = std::get<2>(LDU);
	auto b_new = D.inverse()*b;
	vec r_0 = b - A*x_0;
	double res = 0;
	if (criterion == consecutive)
	{
		vec x_new = x_0;
		vec x_old = x_0-vec::Ones(x_0.size());
		while ((x_new-x_old).norm() > tolerance)
		{
			x_old = x_new;
			x_new = -D.inverse() * (L*x_old + U*x_old) + b_new;
			ic = ic + 1;
		}
		x = x_new;
		std::cout << "norm" << (b - A*x).norm() << std::endl;
	}

	else if (criterion == residual)
	{
		
		vec r = r_0;
		while (r.norm() >= tolerance)
		{
			x = -D.inverse() * (L*x + U*x) + b_new;
			r = b - A * x;
			ic = ic + 1;
		}
	}

	return std::make_tuple(x, ic, (b-A*x).norm());

}



//Gauss_seidel method
std::tuple <vec, int >
gauss_seidel(const mat & A, const vec & b, const vec & x_0,
	const double tolerance,
	const StoppingCriterion criterion)
{
	vec x = x_0;
	auto ic = 0;
	auto LDU = decompose_LDU(A);
	auto L = std::get<0>(LDU);
	auto D = std::get<1>(LDU);
	auto U = std::get<2>(LDU);
	auto b_new = forward_subst(D+L,b);


	if (criterion == consecutive)
	{
		vec x_old = x_0;
		vec diff = x_0;
		diff(0) = 10;
		while ((diff).norm() >= tolerance)
		{
			x = -forward_subst(D+L, U*x_old)+b_new;
			diff = x - x_old;
			x_old = x;
			ic = ic + 1;
		}
		std::cout << std::endl << "norm" << diff.norm() << std::endl;
	}

	else if (criterion == residual)
	{
		vec r_0 = b - A*x_0;
		vec r = r_0;
		while (r.norm() >= tolerance)
		{
			x = -forward_subst(D + L, U*x) + b_new;
			r = b - A * x;
			ic = ic + 1;
		}
	}

	return std::make_tuple(x, ic);
}



std::tuple <vec, int,double>
gauss_seidel(const mat & A, const vec & b,
	const double tolerance,
	const StoppingCriterion criterion)
{
	vec x_0 = vec::Zero(b.size());
	vec x = x_0;
	auto ic = 0;
	auto LDU = decompose_LDU(A);
	auto L = std::get<0>(LDU);
	auto D = std::get<1>(LDU);
	auto U = std::get<2>(LDU);
	auto b_new = forward_subst(D + L, b);
	double res = 0;

	if (criterion == consecutive)
	{
		vec x_old = x_0;
		vec diff = x_0;
		res = 1;
		while (res >= tolerance)
		{
			x = -forward_subst(D + L, U*x_old) + b_new;
			diff = x - x_old;
			x_old = x;
			ic = ic + 1;
			res = diff.norm();
		}
		std::cout << std::endl << "norm" << (b - A * x).norm() << std::endl;
		
	}

	else if (criterion == residual)
	{
		vec r_0 = b - A*x_0;
		vec r = r_0;
		while (r.norm() >= tolerance)
		{
			x = -forward_subst(D + L, U*x) + b_new;
			r = b - A * x;
			ic = ic + 1;
		}
	}

	return std::make_tuple(x, ic, (b - A * x).norm());
}



//sor method
std::tuple <vec, int >
sor(const mat & A, const vec & b, const vec & x_0,
	const double tolerance,
	const StoppingCriterion criterion)
{
	vec x = x_0;
	auto ic = 0;
	auto w = 1.05;
	auto LDU = decompose_LDU(A);
	auto L = std::get<0>(LDU);
	auto D = std::get<1>(LDU);
	auto U = std::get<2>(LDU);

	auto b_new = w * forward_subst(D + w*L, b);


	if (criterion == consecutive)
	{
		vec x_old = x_0;
		vec diff = x_0 - vec::Ones(x_0.size());
		while ((diff).norm() > tolerance)
		{
			x = forward_subst(D + w*L, (1 - w)*D*x_old - w*U*x_old) + b_new;
			diff = x - x_old;
			x_old = x;
			ic = ic + 1;
		}
		std::cout << std::endl << "norm" << diff.norm() << std::endl;
	}

	else if (criterion == residual)
	{
		vec r_0 = b - A*x_0;
		vec r = r_0;
		while (r.norm() >= tolerance)
		{
			x = forward_subst(D + w*L, (1 - w)*D*x - w*U*x) + b_new;
			r = b - A * x;
			ic = ic + 1;
		}
	}

	return std::make_tuple(x, ic);
}


std::tuple <vec, int, double>
sor(const mat & A, const vec & b,
	const double tolerance,
	const StoppingCriterion criterion)
{
	vec x_0 = vec::Zero(b.size());
	vec x = x_0;

	auto ic = 0;
	auto w = 1.05;
	auto LDU = decompose_LDU(A);
	auto L = std::get<0>(LDU);
	auto D = std::get<1>(LDU);
	auto U = std::get<2>(LDU);
	vec b_new = w * forward_subst(D + w*L, b);

	if (criterion == consecutive)
	{
		vec x_old = x_0;
		vec diff = x_0-vec::Ones(x_0.size());
		while ((diff).norm() >= tolerance)
		{
			x = forward_subst(D + w*L, (1 - w)*D*x_old - w*U*x_old) + b_new;
			diff = x - x_old;
			x_old = x;
			ic = ic + 1;
		}
		std::cout << (x - x_old).norm() << std::endl;
	}

	else if (criterion == residual)
	{
		vec r_0 = b - A*x_0;
		vec r = r_0;
		while (r.norm() >= tolerance)
		{
			x = forward_subst(D + w*L, (1 - w)*D*x - w*U*x) + b_new;
			r = b - A * x;
			ic = ic + 1;
		}
	}
	std::cout << std::endl << "norm" << (b - A * x).norm() << std::endl;
	return std::make_tuple(x, ic,(b - A * x).norm());
}