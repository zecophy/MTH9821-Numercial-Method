#pragma once
#include"BinomialTree.h"

ValueAndGreeks AverageBinomialTree(Option& op, int N = 100) {
	auto a = BinomialTree(op, N);
	auto b = BinomialTree(op, N + 1);
	return { (a.value + b.value) / 2,(a.delta + b.delta) / 2,(a.gamma + b.gamma) / 2,
		(a.theta + b.theta) / 2 };
}

