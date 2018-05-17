// Copyright (c) 2018 Denis Skripov aka nizhikebinesi
//
// метод отражений с выбором ведущего столбца
#pragma once

#include <iostream>

#include "numeric_methods.h"
#include "Mygen\Mygen.cpp"

using namespace std;
using namespace numeric_methods;


void get_symmetric(double** a, double** a_inv, int n, double alpha, double beta) {
	mygen(a, a_inv, n, alpha, beta, 1, 2, 0, 1);
}

void get_simple_structure(double** a, double** a_inv, int n, double alpha, double beta) {
	mygen(a, a_inv, n, alpha, beta, 1, 2, 1, 1);
}

void get_jordan_cell(double** a, double** a_inv, int n, double alpha, double beta) {
	mygen(a, a_inv, n, alpha, beta, 0, 0, 2, 1);
}

void test() {
	int cond = 3;

	double** A = create_matrix(N);
	double** A_copy = create_matrix(N);
	double** A_inv = create_matrix(N);
	double** error = create_matrix(N);
	double** residual = create_matrix(N);
	double** E = create_unit_matrix(N);

	// alpha, beta;
	for (double alpha = 0.1; alpha > 1e-17; alpha /= 10) {
		for (double beta = 1.0; beta < 2.0; beta++) {
			switch (cond){
				case 1: {
					get_symmetric(A, A_inv, N, alpha, beta);
					break;
				}
				case 2: {
					get_simple_structure(A, A_inv, N, alpha, beta);
					break;
				}
				/*case 3: {
					get_jordan_cell(A, A_inv, size, alpha, beta);
					break;
				}*/
				default: {
					get_jordan_cell(A, A_inv, N, alpha, beta);
					break;
				}
			}

			cout << "alpha = " << alpha << "; beta = " << beta << '\n';
			
			A_copy = copy_matrix(A, N);
			yet_another_reflection(A, N);
			
			error = matrix_sub(A, A_inv, N);
			residual = matrix_sub(matrix_mult(A_copy, A, N), E, N);
			
			double z, r, zeta;
			z = infinity_norm(error, N);
			r = infinity_norm(residual, N);
			zeta = z / infinity_norm(A_inv, N);

			cout << "||z|| = " << z << "\n";
			cout << "||r|| = " << r << "\n";
			cout << "zeta = " << zeta << "\n";
			cout << "\n\n\n";
		}
	}

	delete_matrix(A, N);
	delete_matrix(A_inv, N);
	delete_matrix(A_copy, N);
	delete_matrix(E, N);
	delete_matrix(error, N);
	delete_matrix(residual, N);
}

int main() {
	/*double** matr = nullptr;
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

	cout << 

	delete_matrix(E, size);
	delete_matrix(copy_matr, size);
	delete_matrix(matr, size);
	*/

	test();

	for (int i = 0; i < 3; i++) {
		wait();
	}



	return 0;
}