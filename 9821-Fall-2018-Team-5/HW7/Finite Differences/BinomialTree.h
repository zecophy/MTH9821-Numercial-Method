#pragma once
#include"Option.h"


double BinomialTree(Option op, int N = 2500) {
	double dt = op.T / N;
	double disc = std::exp(-op.r*dt);
	double u = std::exp(op.sigma*std::sqrt(dt));
	double d = 1 / u;
	double pu = (std::exp((op.r - op.q)*dt) - d) / (u - d);
	double pd = 1 - pu;

	double Smax = op.S0 * std::pow(u, N);
	std::vector<double> V(N + 1);
	for (int i = 0;i <= N;i++) {
		if (op.cp == put)
			V[i] = std::max(op.K - Smax, 0.0);
		else if (op.cp == call)
			V[i] = std::max(-op.K + Smax, 0.0);
		Smax /= (u*u);
	}
	
	// if american option
	
	for (int k = N - 1;k >= 0;k--) {
		double Smax2 = op.S0*std::pow(u, k);
		for (int i = 0;i <= k;i++) {
			// if put
			if (op.cp == put) {
				V[i] = std::max(disc*(pu*V[i] + pd * V[i + 1]), std::max(op.K - Smax2, 0.0));
			}
			//if call
			else if (op.cp == call) {
				V[i] = std::max(disc*(pu*V[i] + pd * V[i + 1]), std::max(-op.K + Smax2, 0.0));
			}
			Smax2 /= (u*u);
		}
	}
	
	return V[0];
}
