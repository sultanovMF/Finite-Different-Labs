#pragma once 

#include <vector>
#include <functional>

struct CauchyProblem
{
	std::function<std::vector<double>(double, std::vector<double>)> foo;
	std::vector<double> y;
	double t0;
	double tf;
};

auto euler_method(const CauchyProblem& problem, double h, int n) {
	auto t = problem.t0;
	auto y = problem.y;

	std::vector<std::pair<double, std::vector<double>>> result;
	result.push_back(std::pair(t, y)); // handler
	for (int i = 1; i < n; i++)
	{
		auto y_prev = problem.foo(t, y);
		for (int j = 0; j < y.size(); ++j)
		{
			y[j] += h * y_prev[j];
		}
		t += h;

		result.push_back(std::pair(t, y)); // handler
	}
	return result;
}

auto adams_method(const CauchyProblem& problem, double h, int n) {
	auto t = problem.t0;
	auto y = problem.y;

	std::vector<std::pair<double, std::vector<double>>> result;
	result.push_back(std::pair(t, y)); // handler T0 Y0

	// Set Y0
	auto y_prev_prev = y;
	auto foo_prev_prev = problem.foo(t, y_prev_prev); // f(T0, Y0)

	// Calculate Y1 using Euler method
	auto y_prev = y_prev_prev;
	for (int j = 0; j < y.size(); ++j)
	{
		y_prev[j] += h * foo_prev_prev[j];
	}

	t += h;
	auto foo_prev = problem.foo(t, y_prev); // f(T1, Y1)

	result.push_back(std::pair(t, y_prev)); // handler T1 Y1

	for (int i = 2; i < n; i++)
	{
		// Calculate Yn+1
		y = y_prev;
		for (int j = 0; j < y.size(); ++j)
		{
			y[j] += h / 2. * (3 * foo_prev[j] - foo_prev_prev[j]);
		}

		t += h;

		result.push_back(std::pair(t, y)); // handler

		y_prev_prev = y_prev;
		y_prev = y;

		foo_prev_prev = foo_prev; // f(Tn, Yn)
		foo_prev = problem.foo(t, y);	  // f(Tn+1, Yn+1)
	}
	return result;
}

auto rk4_method(const CauchyProblem& problem, double h, int n) {
	auto t = problem.t0;
	auto y = problem.y;
	std::vector<std::pair<double, std::vector<double>>> result;

	result.push_back(std::pair(t, y)); // handler T0 Y0

	auto shift = [](double new_h, std::vector<double> k, std::vector<double> y)
	{
		std::vector<double> temp = y;
		for (int i = 0; i < y.size(); ++i)
		{
			temp[i] += new_h / 2. * k[i];
		}
		return temp;
	};

	for (int i = 1; i < n; ++i)
	{
		auto k1 = problem.foo(t, y);
		auto k2 = problem.foo(t + h / 2., shift(h / 2., k1, y));
		auto k3 = problem.foo(t + h / 2., shift(h / 2., k2, y));
		auto k4 = problem.foo(t + h,      shift(h,      k3, y));

		for (int j = 0; j < y.size(); ++j)
		{
			y[j] += h / 6. * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
		}

		t += h;

		result.push_back(std::pair(t, y)); // handler
	}

	return result;
}
