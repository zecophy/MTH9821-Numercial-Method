#pragma once
#include"TBS.h"

ValueAndGreeks TBSR(Option& op, int N = 100) {
	auto a = TBS(op, N);
	auto b = TBS(op, N / 2);
	return { (2 * a.value - b.value),(2 * a.delta - b.delta),(2 * a.gamma - b.gamma),(2 * a.theta - b.theta) };
}
