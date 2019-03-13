#pragma once
#pragma once
#include<cmath>
#include<vector>
#include<algorithm>
#include<iostream>

enum Type { european, american };
enum CP { put, call };

struct ValueAndGreeks {
	double value;
	double delta;
	double gamma;
	double theta;
};

struct Option {
	Type type;
	CP cp;
	double S0;
	double K;
	double T;
	double sigma;
	double r;
	double q;
};

ValueAndGreeks TrinomialTree(Option& op, int N = 100) {
	double dt = op.T / N;
	double disc = std::exp(-op.r*dt);
	double u = std::exp(op.sigma*std::sqrt(3*dt));
	double m = 1;
	double d = 1 / u;
	double pu = 1.0 / 6.0 + (op.r - op.q - 0.5*op.sigma*op.sigma)*std::sqrt(dt / 12 / op.sigma / op.sigma);
	double pm = 2.0 / 3.0;
	double pd = 1 - pu - pm;

	double Smax = op.S0 * std::pow(u, N);
	std::vector<double> V(2 * N + 1);
	for (int i = 0;i <= 2*N;i++) {
		if (op.cp == put)
			V[i] = std::max(op.K - Smax, 0.0);
		else if (op.cp == call)
			V[i] = std::max(-op.K + Smax, 0.0);
		Smax /= u;
	}

	std::vector<std::vector<double>> store{ { 0 },{ 0,0,0 },{ 0,0,0,0,0 } };

	// if european option
	if (op.type == european) {
		for (int k = N - 1;k >= 0;k--) {
			for (int i = 0;i <= 2*k;i++) {
				V[i] = disc * (pu*V[i] + pm * V[i + 1] + pd * V[i + 2]);
				if (k <= 2)store[k][i] = V[i];
			}
		}
	}

	// if american option
	else if (op.type == american) {
		for (int k = N - 1;k >= 0;k--) {
			double Smax2 = op.S0*std::pow(u, k);
			for (int i = 0;i <= 2*k;i++) {
				// if put
				if (op.cp == put) {
					V[i] = std::max(disc * (pu*V[i] + pm * V[i + 1] + pd * V[i + 2]), 
						std::max(op.K - Smax2, 0.0));
				}
				//if call
				else if (op.cp == call) {
					V[i] = std::max(disc * (pu*V[i] + pm * V[i + 1] + pd * V[i + 2]), 
						std::max(-op.K + Smax2, 0.0));
				}
				if (k <= 2)store[k][i] = V[i];
				Smax2 /= u;
			}
		}
	}

	double delta = (store[1][0] - store[1][2]) / (op.S0*u - op.S0*d);
	double gamma = ((store[2][0] - store[2][2]) / (op.S0*u*u - op.S0)
		- (store[2][2] - store[2][4]) / (op.S0 - op.S0*d*d)) / (op.S0*u - op.S0*d);
	double theta = (store[1][1] - store[0][0]) /  dt;
	return { V[0],delta,gamma,theta };
}
