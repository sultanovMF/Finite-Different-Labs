#include <iostream>
#include <span>
#include <vector>
#include <functional>
#include <fstream>
#include <tuple>
#include "matrix_solvers.hpp"

const double EPS = 10e-6;

auto euler_method(std::function<double(double, double)> foo, double y, double t, double h, int n) {
	std::vector<std::pair<double, double>> result;
    for (int i = 0; i < n; i++) {
        y += h * foo(y, t);
        t += h;

		result.push_back(std::pair(t, y)); // handler
    }

	return result;
}

int main() {
	auto foo = [](double y, double t) { return 6*t*t + 5*t*y; };
	auto res = euler_method(foo, 1, 1, 0.01, 10);

	std::fstream output;

	output.open("../../../results/lab_1/ex1.txt", std::ios::out);
	if (!output) {
		std::cout << "File not created!";
	}
	else {
		for (auto [t, y] : res) {
			output << t << "\t" << y << "\n";
		}
		output.close();
	}
    return 0;
}
