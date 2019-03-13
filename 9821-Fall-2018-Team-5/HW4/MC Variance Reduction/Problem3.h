#pragma once
#include"BasketOptionPathDependent.h"

void problem3() {
	Option op1(European, put, 25, 50, 0.5, 0.3, 0, 0.05);
	Option op2(European, put, 30, 50, 0.5, 0.2, 0, 0.05);
	BasketOptionPathDependentPricer pricer(op1, op2, 0.3, 150);

	ofstream out;
	out.open("3_data.csv");

	for (int N = 50;N <= 25600;N *= 2) {
		cout << N << endl;
		pricer.PricerGo(N);
		out << setprecision(10) << N << "," << pricer.value << endl;
	}

	out.close();
	
	return;

}