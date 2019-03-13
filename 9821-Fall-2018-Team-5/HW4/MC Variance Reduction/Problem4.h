#pragma once
#include"HestonOption.h"

void problem4() {
	HestonOption op(European, put, 50, 50, 0.5, 0, 0.05, 0.09, 4, 0.35*0.35, 0.25, -0.15);
	HestonOptionPricer pricer(op, 175);

	ofstream out;
	out.open("4_data.csv");

	for (int N = 500;N <= 16000;N *= 2) {
		cout << N << endl;
		pricer.PricerGo(N);
		out << setprecision(10) << N << "," << pricer.value << "," 
			<< pricer.imlVol << endl;
	}

	out.close();
	return;
}