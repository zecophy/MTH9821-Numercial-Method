//hw8_team5 main.cpp

#include <iostream>
#include <Eigen/Dense>
#include "header.hpp"
#include <fstream>
#include <iomanip>

typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;
<<<<<<< Updated upstream
typedef Eigen::PermutationMatrix <-1, -1, int> permutation;


int main()
{
	std::ofstream out;
	out.open("2_data.csv");
	/*
	//********Test for forward_subst********
	mat fwd(7, 7);
	fwd = mat::Zero(7, 7);
	fwd << -1.73, 0,	0,	0,	0,	0,	0,
		-0.01,0.21,0,0,0,0,0,
		-0.28,	0.76,	0.62,	0,	0,	0,	0,
		1.02,	0.74, - 0.62, - 0.14,	0,	0,	0,
		-0.57,	0.09, - 1.4, - 0.69, - 0.68,	0,	0,
		0.38, - 0.28, - 0.81, - 0.27, - 0.44, - 0.02,	0,
		-1.4, - 1.38, - 0.14, - 2.6,	1.54,	1.65,	0.42;
	vec fwd_b(7);
	fwd_b << 3.63,
		0.42,
		2.15,
		- 1,
		0.46,
		- 1.7,
		- 3.24
		;
	auto fwd_x = forward_subst(fwd, fwd_b);
	std::cout<<fwd_x<<std::endl;
	std::cout<<std::endl;


	out << "fwd" << std::endl;
	out << std::setprecision(10) << fwd_x<<std::endl;
	




	//********Test for backward_subst******************
	mat bwd(5, 5);
	bwd = mat::Zero(5, 5);
	bwd << 1.3, - 4.9,	6.6,	2.4, - 11.3,
		0, - 2.7, - 14.7, - 13.5,	2.2,
		0,	0,	19.1,	0.8, - 9.1,
		0,	0,	0,	0.5,	10.5,
		0,	0,	0,	0,	15.4;
	vec bwd_b(5);
	bwd_b << 1.56,
		- 0.71,
		- 1.1,
		- 1.65,
		- 0.36;
	auto bwd_x = backward_subst(bwd, bwd_b);
	std::cout << bwd_x << std::endl;
	std::cout << std::endl;
	std::cout<<std::endl;

	out << "bwd" << std::endl;
	out << std::setprecision(10) << bwd_x << std::endl;

	//********lu_no_pivoting********
	mat lu_A(6, 6);
	lu_A = mat::Zero(6, 6);
	lu_A << -0.96,	0.82, - 0.88,	0.61,	0.92, - 0.3,
		- 0.76, - 1.17, - 1.25, - 0.56, - 1.82, - 0.05,
		0.73,	0.15,	0.67, - 0.75,	0.1,	0.13,
		- 0.17,	1.72,	0.55, - 1.97, - 2.07,	0.01,
		- 0.39,	1.14, - 0.34, - 1.17,	0.64,	1.06,
		- 0.09, - 0.1, - 0.34,	0.93,	1.6,	0.93
		;
	auto result_n = lu_no_pivoting(lu_A);
	std::cout<<std::endl;
	std::cout<<std::get<0>(result_n) <<std::endl;
	std::cout<<std::endl;
	std::cout<<std::get<1>(result_n) <<std::endl;

	out << "lu_no_pivoting" << std::endl;
	out << std::setprecision(10) << std::get<0>(result_n) << std::endl<<std::endl;
	out << std::setprecision(10) << std::get<1>(result_n) << std::endl;

	//********lu_pivoting********
	mat B(6, 6);
	B << -0.61, - 0.15,- 0.11,	1.99,	1.46, - 1.51,
		0.35, - 0.13, - 1.48,	1.19,	0.1, - 1.4,
		0.68,	1.49, - 1.26,	0.92,	0.17,	0.43,
		- 0.58, - 0.29, - 0.33, - 1.1, - 0.86, - 0.09,
		1.53,	1.27, - 1.18, - 0.06,	2.33,	0.44,
		0.24,	0, - 0.72, - 0.09,	0.41, - 0.49
		;

	auto result_p = lu_row_pivoting(B);

	std::cout<<std::endl<<"********lu_pivoting********"<<std::endl;
	std::cout<<std::get<1>(result_p) <<std::endl;
	std::cout<<std::endl;
	std::cout<<std::get<2>(result_p) <<std::endl;

	mat permu(6, 6);
	permu = mat::Identity(6, 6);
	permu = std::get<0>(result_p) * permu;
	out << "lu_pivoting" << std::endl;
	out << std::setprecision(10) << std::get<1>(result_p) << std::endl << std::endl;
	out << std::setprecision(10) << std::get<2>(result_p) << std::endl << std::endl;
	out << std::setprecision(10) << permu << std::endl << std::endl;


	//*****cholesky test*******
	std::cout <<"*****cholesky test*******"<< std::endl;
	mat a(5, 5);
	a << 2.75,	1.61,	1.51,	1.22,	1.33,
		1.61,	3.93, - 1.75,	2.27,	2.2,
		1.51, - 1.75,	3.32, - 1.1, - 0.59,
		1.22,	2.27, - 1.1,	2.58,	0.99,
		1.33,	2.2, - 0.59,	0.99,	2.57
		;
	vec cho_b(5);
	cho_b << -0.02,
		0.05,
		-0.45,
		-0.19,
		-1.84;

	std::cout << "The matrix A is" << std::endl << a << std::endl;
	auto U1 = cholesky(a);
	std::cout << "The Cholesky factor U is" << std::endl << U1 << std::endl;

	out << "cholesky" << std::endl;
	out << std::setprecision(10) << U1 << std::endl << std::endl;
	out << std::setprecision(10) << cholesky_linear(a,cho_b) << std::endl << std::endl;
	*/

	
	//**********Jacobi****************
	vec m(4);
	m << 16,
		11,
		16,
		19;
	mat N(4, 4);
	N << 5, - 2,	1,	1,
		2,	7,	2, - 1,
		1,	2,	7,	1,
		- 1,	1,	2,	8
		;
	std::cout << std::endl;

	
	std::tuple<vec,int, double> result_j_consecutive = jacobi(N, m, 0.000001, consecutive);

	std::cout << "**********Jacobi****************" << std::endl;
	std::cout << std::get<0>(result_j_consecutive)<<std::endl;
	std::cout << std::get<1>(result_j_consecutive) << std::endl;

	out <<std::setprecision(10) << std::get<0>(result_j_consecutive) << std::endl;
	out << std::setprecision(10) << std::get<1>(result_j_consecutive) << std::endl;
	out << std::setprecision(10) << std::get<2>(result_j_consecutive) << std::endl;

	//**********Gauss_Siedel**********
	vec m_gs(5);
	m_gs << 16,
		11,
		16,
		19,
		20;

	mat N_gs(5, 5);
	N_gs << 10,	3, - 2,	1,	1,
		2,	9,	2,	2, - 1,
		1,	2,	13,	6,	1,
		- 1,	1,	2,	24,	5,
		3,	3,	2, - 4,	19
		;
	std::cout << std::endl;
	std::tuple<vec, int, double> result_gs_consecutive = gauss_seidel(N_gs, m_gs, 0.000001, consecutive);
	

	std::cout << "**********Gauss_Siedel**********" << std::endl;
	std::cout << std::get<0>(result_gs_consecutive) << std::endl;
	std::cout << std::get<1>(result_gs_consecutive) << std::endl;
	out << std::setprecision(10) << std::get<0>(result_gs_consecutive) << std::endl;
	out << std::setprecision(10) << std::get<1>(result_gs_consecutive) << std::endl;
	out << std::setprecision(10) << std::get<2>(result_gs_consecutive) << std::endl;

	//**********SOR*******************
	std::cout << std::endl;

	vec m_sor(6);
	m_sor << 16,
		11,
		16,
		19,
		20,
		-19;

	mat N_sor(6, 6);
	N_sor << 15,	5,	3, - 2,	1,	1,
		2,	26,	6,	2,	2, - 1,
		1,	2,	18,	3,	6,	1,
		- 1,	1,	2,	28,	5,	5,
		3,	3,	2, - 4,	36,	8,
		- 9,	8,	3,	1,	2,	45;


	std::tuple<vec, int, double> result_sor_consecutive = sor(N_sor, m_sor, 0.000001, consecutive);

	std::cout << "**********SOR**********" << std::endl;
	std::cout << std::get<0>(result_sor_consecutive) << std::endl;
	std::cout << std::get<1>(result_sor_consecutive) << std::endl;
	out << std::setprecision(10) << std::get<0>(result_sor_consecutive) << std::endl;
	out << std::setprecision(10) << std::get<1>(result_sor_consecutive) << std::endl;
	out << std::setprecision(10) << std::get<2>(result_sor_consecutive) << std::endl;

	std::cout << std::setprecision(10) << (m_sor - N_sor*std::get<0>(result_sor_consecutive)).norm() << std::endl;
	out.close();
	
	system("pause");
	return 0;
=======
typedef Eigen::PermutationMatrix < -1, -1, int > permutation;


int main() {
    std::ofstream out;
    out.open ("2_data.csv");
    /*
    //********Test for forward_subst********
    mat fwd(7, 7);
    fwd = mat::Zero(7, 7);
    fwd << -1.73, 0,	0,	0,	0,	0,	0,
    	-0.01,0.21,0,0,0,0,0,
    	-0.28,	0.76,	0.62,	0,	0,	0,	0,
    	1.02,	0.74, - 0.62, - 0.14,	0,	0,	0,
    	-0.57,	0.09, - 1.4, - 0.69, - 0.68,	0,	0,
    	0.38, - 0.28, - 0.81, - 0.27, - 0.44, - 0.02,	0,
    	-1.4, - 1.38, - 0.14, - 2.6,	1.54,	1.65,	0.42;
    vec fwd_b(7);
    fwd_b << 3.63,
    	0.42,
    	2.15,
    	- 1,
    	0.46,
    	- 1.7,
    	- 3.24
    	;
    auto fwd_x = forward_subst(fwd, fwd_b);
    std::cout<<fwd_x<<std::endl;
    std::cout<<std::endl;


    out << "fwd" << std::endl;
    out << std::setprecision(10) << fwd_x<<std::endl;





    //********Test for backward_subst******************
    mat bwd(5, 5);
    bwd = mat::Zero(5, 5);
    bwd << 1.3, - 4.9,	6.6,	2.4, - 11.3,
    	0, - 2.7, - 14.7, - 13.5,	2.2,
    	0,	0,	19.1,	0.8, - 9.1,
    	0,	0,	0,	0.5,	10.5,
    	0,	0,	0,	0,	15.4;
    vec bwd_b(5);
    bwd_b << 1.56,
    	- 0.71,
    	- 1.1,
    	- 1.65,
    	- 0.36;
    auto bwd_x = backward_subst(bwd, bwd_b);
    std::cout << bwd_x << std::endl;
    std::cout << std::endl;
    std::cout<<std::endl;

    out << "bwd" << std::endl;
    out << std::setprecision(10) << bwd_x << std::endl;

    //********lu_no_pivoting********
    mat lu_A(6, 6);
    lu_A = mat::Zero(6, 6);
    lu_A << -0.96,	0.82, - 0.88,	0.61,	0.92, - 0.3,
    	- 0.76, - 1.17, - 1.25, - 0.56, - 1.82, - 0.05,
    	0.73,	0.15,	0.67, - 0.75,	0.1,	0.13,
    	- 0.17,	1.72,	0.55, - 1.97, - 2.07,	0.01,
    	- 0.39,	1.14, - 0.34, - 1.17,	0.64,	1.06,
    	- 0.09, - 0.1, - 0.34,	0.93,	1.6,	0.93
    	;
    auto result_n = lu_no_pivoting(lu_A);
    std::cout<<std::endl;
    std::cout<<std::get<0>(result_n) <<std::endl;
    std::cout<<std::endl;
    std::cout<<std::get<1>(result_n) <<std::endl;

    out << "lu_no_pivoting" << std::endl;
    out << std::setprecision(10) << std::get<0>(result_n) << std::endl<<std::endl;
    out << std::setprecision(10) << std::get<1>(result_n) << std::endl;

    //********lu_pivoting********
    mat B(6, 6);
    B << -0.61, - 0.15,- 0.11,	1.99,	1.46, - 1.51,
    	0.35, - 0.13, - 1.48,	1.19,	0.1, - 1.4,
    	0.68,	1.49, - 1.26,	0.92,	0.17,	0.43,
    	- 0.58, - 0.29, - 0.33, - 1.1, - 0.86, - 0.09,
    	1.53,	1.27, - 1.18, - 0.06,	2.33,	0.44,
    	0.24,	0, - 0.72, - 0.09,	0.41, - 0.49
    	;

    auto result_p = lu_row_pivoting(B);

    std::cout<<std::endl<<"********lu_pivoting********"<<std::endl;
    std::cout<<std::get<1>(result_p) <<std::endl;
    std::cout<<std::endl;
    std::cout<<std::get<2>(result_p) <<std::endl;

    mat permu(6, 6);
    permu = mat::Identity(6, 6);
    permu = std::get<0>(result_p) * permu;
    out << "lu_pivoting" << std::endl;
    out << std::setprecision(10) << std::get<1>(result_p) << std::endl << std::endl;
    out << std::setprecision(10) << std::get<2>(result_p) << std::endl << std::endl;
    out << std::setprecision(10) << permu << std::endl << std::endl;


    //*****cholesky test*******
    std::cout <<"*****cholesky test*******"<< std::endl;
    mat a(5, 5);
    a << 2.75,	1.61,	1.51,	1.22,	1.33,
    	1.61,	3.93, - 1.75,	2.27,	2.2,
    	1.51, - 1.75,	3.32, - 1.1, - 0.59,
    	1.22,	2.27, - 1.1,	2.58,	0.99,
    	1.33,	2.2, - 0.59,	0.99,	2.57
    	;
    vec cho_b(5);
    cho_b << -0.02,
    	0.05,
    	-0.45,
    	-0.19,
    	-1.84;

    std::cout << "The matrix A is" << std::endl << a << std::endl;
    auto U1 = cholesky(a);
    std::cout << "The Cholesky factor U is" << std::endl << U1 << std::endl;

    out << "cholesky" << std::endl;
    out << std::setprecision(10) << U1 << std::endl << std::endl;
    out << std::setprecision(10) << cholesky_linear(a,cho_b) << std::endl << std::endl;
    */


    //**********Jacobi****************
    vec m (4);
    m << 16,
    11,
    16,
    19;
    mat N (4, 4);
    N << 5, - 2,	1,	1,
        2,	7,	2, - 1,
        1,	2,	7,	1,
        - 1,	1,	2,	8
        ;
    std::cout << std::endl;


    std::tuple<vec, int, double> result_j_consecutive = jacobi (N, m, 0.000001, consecutive);

    std::cout << "**********Jacobi****************" << std::endl;
    std::cout << std::get<0> (result_j_consecutive) << std::endl;
    std::cout << std::get<1> (result_j_consecutive) << std::endl;

    out << std::setprecision (10) << std::get<0> (result_j_consecutive) << std::endl;
    out << std::setprecision (10) << std::get<1> (result_j_consecutive) << std::endl;
    out << std::setprecision (10) << std::get<2> (result_j_consecutive) << std::endl;

    //**********Gauss_Siedel**********
    vec m_gs (5);
    m_gs << 16,
         11,
         16,
         19,
         20;

    mat N_gs (5, 5);
    N_gs << 10,	3, - 2,	1,	1,
          2,	9,	2,	2, - 1,
          1,	2,	13,	6,	1,
          - 1,	1,	2,	24,	5,
          3,	3,	2, - 4,	19
          ;
    std::cout << std::endl;
    std::tuple<vec, int, double> result_gs_consecutive = gauss_seidel (N_gs, m_gs, 0.000001, consecutive);


    std::cout << "**********Gauss_Siedel**********" << std::endl;
    std::cout << std::get<0> (result_gs_consecutive) << std::endl;
    std::cout << std::get<1> (result_gs_consecutive) << std::endl;
    out << std::setprecision (10) << std::get<0> (result_gs_consecutive) << std::endl;
    out << std::setprecision (10) << std::get<1> (result_gs_consecutive) << std::endl;
    out << std::setprecision (10) << std::get<2> (result_gs_consecutive) << std::endl;

    //**********SOR*******************
    std::cout << std::endl;

    vec m_sor (6);
    m_sor << 16,
          11,
          16,
          19,
          20,
          -19;

    mat N_sor (6, 6);
    N_sor << 15,	5,	3, - 2,	1,	1,
               2,	26,	6,	2,	2, - 1,
               1,	2,	18,	3,	6,	1,
               - 1,	1,	2,	28,	5,	5,
               3,	3,	2, - 4,	36,	8,
               - 9,	8,	3,	1,	2,	45;


    std::tuple<vec, int, double> result_sor_consecutive = sor (N_sor, m_sor, 0.000001, consecutive);

    std::cout << "**********SOR**********" << std::endl;
    std::cout << std::get<0> (result_sor_consecutive) << std::endl;
    std::cout << std::get<1> (result_sor_consecutive) << std::endl;
    out << std::setprecision (10) << std::get<0> (result_sor_consecutive) << std::endl;
    out << std::setprecision (10) << std::get<1> (result_sor_consecutive) << std::endl;
    out << std::setprecision (10) << std::get<2> (result_sor_consecutive) << std::endl;

    std::cout << std::setprecision (10) << (m_sor - N_sor * std::get<0> (result_sor_consecutive) ).norm() << std::endl;
    out.close();
    system ("pause");
    return 0;
>>>>>>> Stashed changes
}


