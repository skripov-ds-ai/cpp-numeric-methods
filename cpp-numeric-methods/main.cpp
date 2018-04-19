// Copyright (c) 2018 Denis Skripov aka nizhikebinesi
//
// метод отражений с выбором ведущего столбца
#pragma once

#include <iostream>

#include "numeric_methods.h"

using namespace std;
using namespace numeric_methods;


int main() {
	double** matr = nullptr;

	size_t size = first_filling(matr);
	
	double** copy_matr = copy_matrix(matr, size);

	inverse_triangle_matrix(matr, size);
	print_matrix(matr, size);

	reflection_method(copy_matr, size);
	print_matrix(copy_matr, size);

	cout << "\n\n";

	delete_matrix(copy_matr, size);
	delete_matrix(matr, size);

	for (int i = 0; i < 3; i++) {
		wait();
	}

	return 0;
}