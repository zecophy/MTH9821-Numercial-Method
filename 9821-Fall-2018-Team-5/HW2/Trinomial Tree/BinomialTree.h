#pragma once
#include"TrinomialTree.h"

ValueAndGreeks BinomialTree(Option& op, int N = 100) {
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

	std::vector<std::vector<double>> store{ {0},{0,0},{0,0,0} };

	// if european option
	if (op.type == european){
		for (int k = N - 1;k >= 0;k--) {
			for (int i = 0;i <= k;i++) {
				V[i] = disc * (pu*V[i] + pd * V[i + 1]);
				if (k <= 2)store[k][i] = V[i];
			}
		}
	}

	// if american option
	else if (op.type == american) {
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
				if (k <= 2)store[k][i] = V[i];
				Smax2 /= (u*u);
			}
		}
	}
	
	double delta = (store[1][0] - store[1][1]) / (op.S0*u - op.S0*d);
	double gamma = ((store[2][0] - store[2][1]) / (op.S0*u*u - op.S0)
		- (store[2][1] - store[2][2]) / (op.S0 - op.S0*d*d)) / ((op.S0*u*u - op.S0*d*d) / 2);
	double theta = (store[2][1] - store[0][0]) / 2 / dt;
	return { V[0],delta,gamma,theta };
}
