// Copyright (c) 2018 Denis Skripov aka nizhikebinesi
//
// метод отражений с выбором ведущего столбца
#include <stdio.h>
#include <math.h>
#include <cstdlib>

#include "numeric_methods.cpp"

using namespace std;
using namespace numeric_methods;

int main() {
	double** matr = nullptr;

	size_t size = first_filling(matr);

	inverse_triangle_matrix(matr, size);

	print_matrix(matr, size);

	delete_matrix(matr, size);

	for (int i = 0; i < 3; i++) {
		wait();
	}

	return 0;
}