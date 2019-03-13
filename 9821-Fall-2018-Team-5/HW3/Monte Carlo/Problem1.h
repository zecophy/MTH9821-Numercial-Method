#pragma once
#include"NonPathDependent.h"

void T1() {
	Option p(European, put, 43, 45, 0.75, 0.3, 0.01, 0.03);
	Option c(European, call, 43, 45, 0.75, 0.3, 0.01, 0.03);
	//p.calculateBS();
	//c.calculateBS();
	int N = 10000;
	shared_ptr<Generator> gen_ptr1 = make_shared<InverseTransform>();
	shared_ptr<Generator> gen_ptr2 = make_shared<InverseTransform>();
	NonPathDependent pricer_P(N, gen_ptr1, p);
	NonPathDependent pricer_C(N, gen_ptr2, c);
	vector<NonPathDependent> pricers{ pricer_C,pricer_P };
	pricer_C.generateVegas();


	std::ofstream out;
	out.open("1_data.csv");

	for (auto pricer : pricers) {
		dynamic_pointer_cast<InverseTransform,Generator>(pricer.gen)->U.seed=1;
		N = 10000;
		vector<double> line1, line2, line3, line4, line5, line6, line7;
		for (;N <= 5120000;N *= 2) {
			cout << N << endl;
			line1.push_back(N);
			pricer.N = N;

			auto temp1 = pricer.MCvalue();
			line2.push_back(temp1);
			line3.push_back(abs(temp1 - pricer.op.BSvalue)*sqrt(N));

			auto temp2 = pricer.MCdelta();
			line4.push_back(temp2);
			line5.push_back(abs(temp2 - pricer.op.BSdelta)*sqrt(N));

			auto temp3 = pricer.MCvega();
			line6.push_back(temp3);
			line7.push_back(abs(temp3 - pricer.op.BSvega)*sqrt(N));
		}

		out << setprecision(10) << "BS value," << pricer.op.BSvalue
			<< ",BS delta," << pricer.op.BSdelta
			<< ",BS vega," << pricer.op.BSvega << endl;

		vector<vector<double>> part1 = { line1,line2,line3,line4,line5,line6,line7 };
		for (int i = 0;i < 7;i++) {
			out << setprecision(10);
			auto temp = part1[i];
			for (auto item : temp) {
				out << item << ",";
			}
			out << endl;
		}
		out << endl << endl;
	}

	//system("pause");
	return;
}
