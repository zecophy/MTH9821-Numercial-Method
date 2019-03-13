#include"TBSR.h"
#include"AverageBinomialTree.h"
#include <fstream>
#include<iomanip>

int main() {
	Option eu{ european,put,41,39,1,0.285,0.03,0.005 };
	Option am{ american,put,41,39,1,0.285,0.03,0.005 };
	/*
	
	// for the first problem
	std::ofstream out;
	out.open("1_data.csv");
	auto bs = BS(eu);
	out <<std::setprecision(10)<< "BS value," << bs.value << ",BS delta," << bs.delta << ",BS gamma," << bs.gamma
		<< ",BS theta," << bs.theta << std::endl;
	out << std::endl;
	for (int N = 10;N <= 1280;N *= 2) {
		std::cout << "N: " << N << std::endl;
		auto temp = TrinomialTree(eu, N);
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
		auto temp = TBS(eu, N);
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
		auto temp = TBSR(eu, N);
		out <<std::setprecision(10)
			<< N << "," << temp.value << "," << std::abs(temp.value - bs.value) << ","
			<< N * std::abs(temp.value - bs.value) << "," << N * N * std::abs(temp.value - bs.value) << ","
			<< temp.delta << "," << std::abs(temp.delta - bs.delta) << ","
			<< temp.gamma << "," << std::abs(temp.gamma - bs.gamma) << ","
			<< temp.theta << "," << std::abs(temp.theta - bs.theta) << std::endl;
	}
	*/

	//////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////
	
	/*
	// for the second problem
	std::ofstream out;
	out.open("2_data.csv");
	auto bs = AverageBinomialTree(am, 10000);

	out <<std::setprecision(10)<< "BS value," << bs.value << ",BS delta," << bs.delta << ",BS gamma," << bs.gamma
		<< ",BS theta," << bs.theta << std::endl;
	out << std::endl;
	
	for (int N = 10;N <= 1280;N *= 2) {
		std::cout << "N: " << N << std::endl;
		auto temp = TrinomialTree(am, N);
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
		auto temp = TBS(am, N);
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
		auto temp = TBSR(am, N);
		out <<std::setprecision(10)
			<< N << "," << temp.value << "," << std::abs(temp.value - bs.value) << ","
			<< N * std::abs(temp.value - bs.value) << "," << N * N * std::abs(temp.value - bs.value) << ","
			<< temp.delta << "," << std::abs(temp.delta - bs.delta) << ","
			<< temp.gamma << "," << std::abs(temp.gamma - bs.gamma) << ","
			<< temp.theta << "," << std::abs(temp.theta - bs.theta) << std::endl;
	}
	*/

	//////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////
	
	// for the second problem with variance reduction
	
	std::ofstream out;
	out.open("3_data.csv");
	auto bs = AverageBinomialTree(am, 10000);
	out <<std::setprecision(10)<< "BS value," << bs.value << ",BS delta," << bs.delta << ",BS gamma," << bs.gamma
		<< ",BS theta," << bs.theta << std::endl;
	out << std::endl;

	auto bs_eu = BS(eu);


	for (int N = 10;N <= 1280;N *= 2) {
		std::cout << "N: " << N << std::endl;
		auto temp = TrinomialTree(am, N);
		auto temp_eu = TrinomialTree(eu, N);

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
		auto temp = TBS(am, N);
		auto temp_eu = TBS(eu, N);

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
		auto temp = TBSR(am, N);
		auto temp_eu = TBSR(eu, N);

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
			<< std::abs(temp.theta + bs_eu.theta - temp_eu.theta- bs.theta) << std::endl;
	}
	
	return 0;
}