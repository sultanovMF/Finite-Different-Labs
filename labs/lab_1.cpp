#include <iostream>
#include <span>
#include <vector>
#include <functional>
#include <fstream>
#include <tuple>
#include <format> // works in MCVS 17.1.1 and above, for linux wait GCC 13 (͡ ° ͜ʖ ͡ °)
#include <array>
#include "ode.hpp"
#include "matrix_solvers.hpp"
#include <span>

// http://www.ees.nmt.edu/outside/courses/hyd510/PDFs/Lecture%20notes/Lectures%20Part%202.6%20FDMs.pdf

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
			output << t << "\t";
			for (int i = 0; i < y.size(); ++i) {
				output << y[i];
				if (i != y.size() - 1) {
					output << "\t";
				}
			}
			output << "\n";
		}
		output.close();
	}
}

void write_results(std::fstream &output, class std::vector<std::pair<double, double>> res, int exnum, int n) {
	std::string path = std::format("../../../results/lab_1/ex_{}_{}.txt", exnum, n);
	output.open(path, std::ios::out);
	if (!output)
	{
		std::cout << "File not created!";
	}
	else
	{
		for (auto [t, y] : res) {
			output << t << "\t" << y << "\n";
		}
		output.close();
	}
}


struct BVP1d {
	// Second Order BVP problem
	// u'' + p(x) u' + q(x) u = f(x)
	// with mixed boundary conditions

	std::function<double (double)> p;
	std::function<double (double)> q;
	std::function<double (double)> f;
	
	// left boundary condition coeffs.
	double alpha[2];
	double A;

	// right boundary condition coeffs.
	double beta[2];
	double B;

	double x0;
	double xl;
};

auto fd1d (const BVP1d& problem, int n) {
	double h = (problem.xl - problem.x0) / (n - 1.);

	std::vector<std::pair<double, double>> result (n); // TODO тут пара вообще не оправдана, но я экономлю время

	std::vector<double> A(n);
	std::vector<double> B(n);
	std::vector<double> C(n);
	std::vector<double> D(n);

	// TODO как аппрроксимировать-то вторым порядком точности???!?!!
	//left points approximation x0 
	// A[0] = 2 - h * problem.p(xi); // (2 - p(Xi) * h) Ui-1
	// B[0] = 2 * h * h * problem.q(xi) - 4; // (2 q(Xi) h^2 - 4) Ui
	// C[0] = 2 + h * problem.p(xi); //(2 + h p(Xi) ) Ui+1
	// D[0] = 2 * h * h * problem.f(xi); // 2 h^2 f(Xi)

	for (int i = 1; i < n; ++i) {
		// inner points approximation
		double xi = problem.x0 + h * i;

		A[i] = 2 - h *     problem.p(xi);     // (2 - p(Xi) * h) Ui-1
		B[i] = 2 * h * h * problem.q(xi) - 4; // (2 q(Xi) h^2 - 4) Ui
		C[i] = 2 + h *     problem.p(xi);     //(2 + h p(Xi) ) Ui+1
		D[i] = 2 * h * h * problem.f(xi);     // 2 h^2 f(Xi)

		result.push_back(std::pair(xi, 0));
	}
	// right points approximation xl

	// sA[n-1] = 2 - h * problem.p(xi); // (2 - p(Xi) * h) Ui-1
	// sB[n-1] = 2 * h * h * problem.q(xi) - 4; // (2 q(Xi) h^2 - 4) Ui
	// sC[n-1] = 2 + h * problem.p(xi); //(2 + h p(Xi) ) Ui+1
	// sD[n-1] = 2 * h * h * problem.f(xi); // 2 h^2 f(Xi)

	auto y = solve_tridiagonal(A, B, C, D);

	for (int i = 0; i < n; ++i) { // 
		result[i].second = y[i];
	}

	return result;

}

int main()
{
	std::fstream output;
	{
		// Part I
		CauchyProblem problem;
		problem.foo = [] (double t, std::vector<double> y) {
			auto x_t = -y[1]; // dx/dt = v
			auto v_t = y[0];
			return std::vector<double>({x_t, v_t});
		};

		problem.t0 = 0;
		problem.tf = 10;
		problem.y = {1, 1};

		double t_len = (problem.tf - problem.t0);

		for (int n : {10, 100, 1000, 10000})
		{
			double h = t_len / n;
			write_results(output, euler_method(problem, h, n), 1, n);
			write_results(output, adams_method(problem, h, n), 2, n);
			write_results(output, rk4_method  (problem, h, n), 3, n);
		}
	}
	{
		// Part II
		BVP1d problem;
		problem.p = [] (double x) {
			return 0.;
		};

		problem.q = [] (double x) {
			return 0.;
		};
		
		problem.f = [] (double x) {
			return 0.;
		};

		problem.A = 0;
		problem.alpha[0] = 0;
		problem.alpha[1] = 0;

		problem.B = 0;
		problem.beta[0] = 0;
		problem.beta[1] = 0;

		problem.x0 = 0;
		problem.xl = 1;

		for (auto n : {10, 100, 1000, 10000}) {
			write_results(output, fd1d(problem, n), 4, n);
		}
		


	}
	return 0;
}
