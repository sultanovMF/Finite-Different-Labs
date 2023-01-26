#include <iostream>
#include <span>
#include <vector>

#include "matrix_solvers.hpp"

int main(int, char **)
{
	double A[4] = {0, 1, 1, 1};
	double B[4] = {2, 10, -5, 4};
	double C[4] = {1, -5, 2, 0};
	double D[4] = {-5, -18, -40, -27};
	auto X = solve_tridiagonal(A, B, C, D);

	double target_x[] = {-3, 1, 5, -8};
	for (int i = 0; i < 4; ++i)
	{
		std::cout << X[i] << " ";
	}

    return 0;
}
