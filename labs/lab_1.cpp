#include <iostream>
#include <span>
#include <vector>
#include <functional>
#include <fstream>
#include <tuple>
#include <format> // works in MCVS 17.1.1, for linux wait GCC 13 (͡ ° ͜ʖ ͡ °)
#include <array>
#include "matrix_solvers.hpp"

const double EPS = 10e-6;

struct ODEProblem
{
	std::function<std::vector<double>(double, std::vector<double>)> foo;
	std::vector<double> y;
	double t0;
	double tf;
};

auto euler_method(const ODEProblem& problem, double h, int n) {
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

auto adams_method(ODEProblem problem, double h, int n) {
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

auto rk4_method(ODEProblem problem, double h, int n) {
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

void write_results(std::fstream &output, class std::vector<std::pair<double, std::vector<double>>> res, int exnum, int n) {
	std::string path = std::format("../../../results/lab_1/ex_{}_{}.txt", exnum, n);
	output.open(path, std::ios::out);
	if (!output)
	{
		std::cout << "File not created!";
	}
	else
	{
		for (auto [t, y] : res)
		{
			output << t << "\t" << y[0] << "\t" << y[1] << "\n";
		}
		output.close();
	}
}

int main()
{
	{
		ODEProblem problem;
		problem.foo = [] (double t, std::vector<double> y) {
			auto x_t = -y[1]; // dx/dt = v
			auto v_t = y[0];
			return std::vector<double>({x_t, v_t});
		};

		problem.t0 = 0;
		problem.tf = 10;
		problem.y = {1, 1};

		double t_len = (problem.tf - problem.t0);

		std::fstream output;
		for (int n : {10, 100, 1000, 10000})
		{
			double h = t_len / n;
			write_results(output, euler_method(problem, h, n), 1, n);
			write_results(output, adams_method(problem, h, n), 2, n);
			write_results(output, rk4_method  (problem, h, n), 3, n);
		}
	}

	return 0;
}
