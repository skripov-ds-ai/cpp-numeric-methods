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
	size_t size = hand_filling(matr);

	double** copy_matr = copy_matrix(matr, size);

	print_matrix(matr, size);
	cout << "\n\n";

	yet_another_reflection(matr, size);
	print_matrix(matr, size);
	cout << "\n\n";
	double** E = matrix_mult(matr, copy_matr, size);
	cout << "E:\n\n";
	print_matrix(E, size);

	delete_matrix(E, size);
	delete_matrix(copy_matr, size);
	delete_matrix(matr, size);

	for (int i = 0; i < 3; i++) {
		wait();
	}

	return 0;
}