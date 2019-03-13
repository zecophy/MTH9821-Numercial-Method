#pragma once
#include"BBS.h"

ValueAndGreeks BBSR(Option& op, int N = 100) {
	auto a = BBS(op, N);
	auto b = BBS(op, N / 2);
	return { (2 * a.value - b.value),(2 * a.delta - b.delta),(2 * a.gamma - b.gamma),(2 * a.theta - b.theta) };
}