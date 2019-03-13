#include"BBSR.h"
#include <fstream>
#include<iomanip>


int main() {
	Option eu{ european,put,45,41,1,0.285,0.025,0.005 };
	Option am{ american,put,45,41,1,0.285,0.025,0.005 };
	
	// for the first problem
	/*
	std::ofstream out;
	out.open("1_data.csv");
	out << "European,,"<<BS(eu).value<<",,,,American"<<",,"<<BS(am).value<<std::endl;

	out <<std::setprecision(10)<< "N" << "," << "BiTree" << "," << "AvgBiTree" << ","
		<< "BBS" << "," << "BBSR"<<",,";
	out <<std::setprecision(10)<< "BiTree" << "," << "AvgBiTree" << "," << "BBS" << "," << "BBSR" << std::endl;

	for (int N = 10;N <= 100;N += 1) {
		std::cout << "N: " << N << std::endl;
		out <<std::setprecision(10)<< N << "," << BinomialTree(eu, N).value << "," << AverageBinomialTree(eu, N).value << ","
			<< BBS(eu, N).value << "," << BBSR(eu, N).value <<",,";
		out <<std::setprecision(10)<< BinomialTree(am, N).value << "," << AverageBinomialTree(am, N).value << ","
			<< BBS(am, N).value << "," << BBSR(am, N).value << std::endl;
	}
	out.close();
	*/

	////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////
	/*
	
	// for the second problem
	std::ofstream out;
	out.open("2_data.csv");
	auto bs = BS(eu);
	out <<std::setprecision(10)<< "BS value," << bs.value << ",BS delta," << bs.delta << ",BS gamma," << bs.gamma
		<< ",BS theta," << bs.theta << std::endl;
	out << std::endl;
	for (int N = 10;N <= 1280;N *= 2) {
		std::cout << "N: " << N << std::endl;
		auto temp = BinomialTree(eu, N);
		out <<std::setprecision(10)
			<< N << "," << temp.value << "," << std::abs(temp.value - bs.value) << ","
			<< N*std::abs(temp.value-bs.value) << "," << N * N * std::abs(temp.value - bs.value) << ","
			<< temp.delta << "," << std::abs(temp.delta - bs.delta) << ","
			<< temp.gamma << "," << std::abs(temp.gamma - bs.gamma) << ","
			<< temp.theta << "," << std::abs(temp.theta - bs.theta) << std::endl;
	}
	out << std::endl;
	for (int N = 10;N <= 1280;N *= 2) {
		std::cout << "N: " << N << std::endl;
		auto temp = AverageBinomialTree(eu, N);
		out <<std::setprecision(10)
			<< N << "," << temp.value << "," << std::abs(temp.value - bs.value) << ","
			<< N * std::abs(temp.value - bs.value) << "," << N * N * std::abs(temp.value - bs.value) << ","
			<< temp.delta << "," << std::abs(temp.delta - bs.delta) << ","
			<< temp.gamma << "," << std::abs(temp.gamma - bs.gamma) << ","
			<< temp.theta << "," << std::abs(temp.theta - bs.theta) << std::endl;
	}
	out << std::endl;
	for (int N = 10;N <= 1280;N *= 2) {
		std::cout << "N: " << N << std::endl;
		auto temp = BBS(eu, N);
		out <<std::setprecision(10)
			<< N << "," << temp.value << "," << std::abs(temp.value - bs.value) << ","
			<< N * std::abs(temp.value - bs.value) << "," << N * N * std::abs(temp.value - bs.value) << ","
			<< temp.delta << "," << std::abs(temp.delta - bs.delta) << ","
			<< temp.gamma << "," << std::abs(temp.gamma - bs.gamma) << ","
			<< temp.theta << "," << std::abs(temp.theta - bs.theta) << std::endl;
	}
	out << std::endl;
	for (int N = 10;N <= 1280;N *= 2) {
		std::cout << "N: " << N << std::endl;
		auto temp = BBSR(eu, N);
		out <<std::setprecision(10)
			<< N << "," << temp.value << "," << std::abs(temp.value - bs.value) << ","
			<< N * std::abs(temp.value - bs.value) << "," << N * N * std::abs(temp.value - bs.value) << ","
			<< temp.delta << "," << std::abs(temp.delta - bs.delta) << ","
			<< temp.gamma << "," << std::abs(temp.gamma - bs.gamma) << ","
			<< temp.theta << "," << std::abs(temp.theta - bs.theta) << std::endl;
	}

	*/

	///////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////
	/*
	
	// for the third problem
	std::ofstream out;
	out.open("3_data.csv");
	auto bs = AverageBinomialTree(am, 10000);
	out <<std::setprecision(10)<< "BS value," << bs.value << ",BS delta," << bs.delta << ",BS gamma," << bs.gamma
		<< ",BS theta," << bs.theta << std::endl;
	out << std::endl;
	for (int N = 10;N <= 1280;N *= 2) {
		std::cout << "N: " << N << std::endl;
		auto temp = BinomialTree(am, N);
		out <<std::setprecision(10)
			<< N << "," << temp.value << "," << std::abs(temp.value - bs.value) << ","
			<< N * std::abs(temp.value - bs.value) << "," << N * N * std::abs(temp.value - bs.value) << ","
			<< temp.delta << "," << std::abs(temp.delta - bs.delta) << ","
			<< temp.gamma << "," << std::abs(temp.gamma - bs.gamma) << ","
			<< temp.theta << "," << std::abs(temp.theta - bs.theta) << std::endl;
	}
	out << std::endl;
	for (int N = 10;N <= 1280;N *= 2) {
		std::cout << "N: " << N << std::endl;
		auto temp = AverageBinomialTree(am, N);
		out <<std::setprecision(10)
			<< N << "," << temp.value << "," << std::abs(temp.value - bs.value) << ","
			<< N * std::abs(temp.value - bs.value) << "," << N * N * std::abs(temp.value - bs.value) << ","
			<< temp.delta << "," << std::abs(temp.delta - bs.delta) << ","
			<< temp.gamma << "," << std::abs(temp.gamma - bs.gamma) << ","
			<< temp.theta << "," << std::abs(temp.theta - bs.theta) << std::endl;
	}
	out << std::endl;
	for (int N = 10;N <= 1280;N *= 2) {
		std::cout << "N: " << N << std::endl;
		auto temp = BBS(am, N);
		out <<std::setprecision(10)
			<< N << "," << temp.value << "," << std::abs(temp.value - bs.value) << ","
			<< N * std::abs(temp.value - bs.value) << "," << N * N * std::abs(temp.value - bs.value) << ","
			<< temp.delta << "," << std::abs(temp.delta - bs.delta) << ","
			<< temp.gamma << "," << std::abs(temp.gamma - bs.gamma) << ","
			<< temp.theta << "," << std::abs(temp.theta - bs.theta) << std::endl;
	}
	out << std::endl;
	for (int N = 10;N <= 1280;N *= 2) {
		std::cout << "N: " << N << std::endl;
		auto temp = BBSR(am, N);
		out <<std::setprecision(10)
			<< N << "," << temp.value << "," << std::abs(temp.value - bs.value) << ","
			<< N * std::abs(temp.value - bs.value) << "," << N * N * std::abs(temp.value - bs.value) << ","
			<< temp.delta << "," << std::abs(temp.delta - bs.delta) << ","
			<< temp.gamma << "," << std::abs(temp.gamma - bs.gamma) << ","
			<< temp.theta << "," << std::abs(temp.theta - bs.theta) << std::endl;
	}
	*/

	///////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////
	/*
	// for the third problem with variance reduction
	std::ofstream out;
	out.open("4_data.csv");
	auto bs = AverageBinomialTree(am, 10000);
	out <<std::setprecision(10)<< "BS value," << bs.value << ",BS delta," << bs.delta << ",BS gamma," << bs.gamma
	<< ",BS theta," << bs.theta << std::endl;
	out << std::endl;

	auto bs_eu = BS(eu);


	for (int N = 10;N <= 1280;N *= 2) {
		std::cout << "N: " << N << std::endl;
		auto temp = BinomialTree(am, N);
		auto temp_eu = BinomialTree(eu, N);

		out << std::setprecision(10)
			<< N << "," << temp.value + bs_eu.value - temp_eu.value << ","
			<< std::abs(temp.value + bs_eu.value - temp_eu.value - bs.value) << ","
			<< N * std::abs(temp.value + bs_eu.value - temp_eu.value - bs.value) << ","
			<< N * N * std::abs(temp.value + bs_eu.value - temp_eu.value - bs.value) << ","
			<< temp.delta + bs_eu.delta - temp_eu.delta << ","
			<< std::abs(temp.delta + bs_eu.delta - temp_eu.delta - bs.delta) << ","
			<< temp.gamma + bs_eu.gamma - temp_eu.gamma << ","
			<< std::abs(temp.gamma + bs_eu.gamma - temp_eu.gamma - bs.gamma) << ","
			<< temp.theta + bs_eu.theta - temp_eu.theta << ","
			<< std::abs(temp.theta + bs_eu.theta - temp_eu.theta - bs.theta) << std::endl;
	}
	out << std::endl;
	for (int N = 10;N <= 1280;N *= 2) {
		std::cout << "N: " << N << std::endl;
		auto temp = AverageBinomialTree(am, N);
		auto temp_eu = AverageBinomialTree(eu, N);

		out << std::setprecision(10) 
			<< N << "," << temp.value + bs_eu.value - temp_eu.value << ","
			<< std::abs(temp.value + bs_eu.value - temp_eu.value - bs.value) << ","
			<< N * std::abs(temp.value + bs_eu.value - temp_eu.value - bs.value) << ","
			<< N * N * std::abs(temp.value + bs_eu.value - temp_eu.value - bs.value) << ","
			<< temp.delta + bs_eu.delta - temp_eu.delta << ","
			<< std::abs(temp.delta + bs_eu.delta - temp_eu.delta - bs.delta) << ","
			<< temp.gamma + bs_eu.gamma - temp_eu.gamma << ","
			<< std::abs(temp.gamma + bs_eu.gamma - temp_eu.gamma - bs.gamma) << ","
			<< temp.theta + bs_eu.theta - temp_eu.theta << ","
			<< std::abs(temp.theta + bs_eu.theta - temp_eu.theta - bs.theta) << std::endl;
	}
	out << std::endl;
	for (int N = 10;N <= 1280;N *= 2) {
		std::cout << "N: " << N << std::endl;
		auto temp = BBS(am, N);
		auto temp_eu = BBS(eu, N);

		out << std::setprecision(10) 
			<< N << "," << temp.value + bs_eu.value - temp_eu.value << ","
			<< std::abs(temp.value + bs_eu.value - temp_eu.value - bs.value) << ","
			<< N * std::abs(temp.value + bs_eu.value - temp_eu.value - bs.value) << ","
			<< N * N * std::abs(temp.value + bs_eu.value - temp_eu.value - bs.value) << ","
			<< temp.delta + bs_eu.delta - temp_eu.delta << ","
			<< std::abs(temp.delta + bs_eu.delta - temp_eu.delta - bs.delta) << ","
			<< temp.gamma + bs_eu.gamma - temp_eu.gamma << ","
			<< std::abs(temp.gamma + bs_eu.gamma - temp_eu.gamma - bs.gamma) << ","
			<< temp.theta + bs_eu.theta - temp_eu.theta << ","
			<< std::abs(temp.theta + bs_eu.theta - temp_eu.theta - bs.theta) << std::endl;
	}
	out << std::endl;
	for (int N = 10;N <= 1280;N *= 2) {
		std::cout << "N: " << N << std::endl;
		auto temp = BBSR(am, N);
		auto temp_eu = BBSR(eu, N);

		out << std::setprecision(10) 
			<< N << "," << temp.value + bs_eu.value - temp_eu.value << ","
			<< std::abs(temp.value + bs_eu.value - temp_eu.value - bs.value) << ","
			<< N * std::abs(temp.value + bs_eu.value - temp_eu.value - bs.value) << ","
			<< N * N * std::abs(temp.value + bs_eu.value - temp_eu.value - bs.value) << ","
			<< temp.delta + bs_eu.delta - temp_eu.delta << ","
			<< std::abs(temp.delta + bs_eu.delta - temp_eu.delta - bs.delta) << ","
			<< temp.gamma + bs_eu.gamma - temp_eu.gamma << ","
			<< std::abs(temp.gamma + bs_eu.gamma - temp_eu.gamma - bs.gamma) << ","
			<< temp.theta + bs_eu.theta - temp_eu.theta << ","
			<< std::abs(temp.theta + bs_eu.theta - temp_eu.theta - bs.theta) << std::endl;
	}
	*/

	return 0;
}