#pragma once
#include"BinomialTree.h"

void problem3() {
	int n = 0;
	double sigma_0 = 0.1, sigma_1 = 0.5;
	double v_m = 3.90;

	double value_0, value_1;
	value_0 = BinomialTree(Option(American, put, 42, 45, 0.75, sigma_0, 0.02, 0.04));
	value_1 = BinomialTree(Option(American, put, 42, 45, 0.75, sigma_1, 0.02, 0.04));

	double sigma_2;
	while (abs(sigma_1 - sigma_0) > 0.00001) {
		sigma_2 = sigma_1 - (value_1 - v_m) / ((value_1 - value_0) / (sigma_1 - sigma_0));
		sigma_0 = sigma_1;
		sigma_1 = sigma_2;
		value_0 = value_1;
		value_1= BinomialTree(Option(American, put, 42, 45, 0.75, sigma_1, 0.02, 0.04));
		n += 1;
		cout << n << ": " << sigma_2 << endl;
	}

	ofstream out;
	out.open("3_data.csv");
	
	out << setprecision(10) << (sigma_0 + sigma_1)*0.5 << endl;
	out << n << endl;

	out.close();

	return;
}
