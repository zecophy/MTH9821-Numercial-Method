#pragma once
#include"BasketOption.h"

void problem2() {
	Option op1(European, put, 25, 50, 0.5, 0.3, 0, 0.05);
	Option op2(European, put, 30, 50, 0.5, 0.2, 0, 0.05);
	BasketOptionPricer pricer(op1,op2,0.3);

	ofstream out;
	out.open("2_data.csv");

	for (int N = 10000;N <= 2560000;N *= 2) {
		cout << N << endl;
		pricer.PricerGo(N);
		out << setprecision(10) << N << "," << pricer.value << endl;
	}

	out.close();
	return;
}
