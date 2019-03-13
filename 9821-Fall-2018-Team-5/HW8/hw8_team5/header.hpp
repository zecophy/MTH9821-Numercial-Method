#ifndef header_hpp
#define header_hpp

#include <iostream>
#include <Eigen/Dense>

typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;
typedef Eigen::PermutationMatrix <-1, -1, int> permutation;

vec forward_subst(const mat & L, const vec & b);
vec backward_subst(const mat & U, const vec & b);

std::tuple <mat, mat > lu_no_pivoting(mat A);
std::tuple <permutation, mat, mat > lu_row_pivoting(mat A);

mat cholesky(mat A);
vec cholesky_linear(mat &A, vec &b);

enum StoppingCriterion { consecutive, residual };
std::tuple<mat, mat, mat> decompose_LDU(const mat&A);
std::tuple <vec, int,double>
jacobi(const mat & A, const vec & b, const vec & x_0,
	const double tolerance,
	const StoppingCriterion criterion);

std::tuple <vec, int, double >
jacobi(const mat & A, const vec & b,
	const double tolerance,
	const StoppingCriterion criterion);

std::tuple <vec, int >
gauss_seidel(const mat & A, const vec & b, const vec & x_0,
	const double tolerance,
	const StoppingCriterion criterion);

std::tuple <vec, int, double >
gauss_seidel(const mat & A, const vec & b,
	const double tolerance,
	const StoppingCriterion criterion);


std::tuple <vec, int >
sor(const mat & A, const vec & b, const vec & x_0,
	const double tolerance,
	const StoppingCriterion criterion);

std::tuple <vec, int, double >
sor(const mat & A, const vec & b,
	const double tolerance,
	const StoppingCriterion criterion);

#endif