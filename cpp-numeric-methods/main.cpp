// Copyright (c) 2018 Denis Skripov aka nizhikebinesi
//
// метод отражений с выбором ведущего столбца
#pragma once

#include <iostream>

#include "numeric_methods.h"
#include "Mygen\Mygen.cpp"

#include <ctime>

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

void test_2() {
	int cond = 1;
	int err = 0;

	cout << "start 1!\n\n";

	for (int n = 100; n < 120; n += 20) {
		for (double eps = 1e-9; eps < 1 + 1e-9; eps++) {
			int its = 0;

			long start = clock();

			double* exact_sol = gen_exact_phi(n);
			double* f = gen_exact_f(n);

			double* sol = min_error_method_3_17(n, f, err, its, eps);

			double* error = vect_sub(sol, exact_sol, size_3_17(n));

			double* Ax = create_vector(size_3_17(n));
			mult_3_14(n, sol, Ax);

			double* residual = vect_sub(Ax, f, size_3_17(n));

			double z, r, zeta, ro;
			z = vect_norm(error, size_3_17(n));
			r = vect_norm(residual, size_3_17(n));
			zeta = z / vect_norm(exact_sol, size_3_17(n));
			ro = r / vect_norm(f, size_3_17(n));

			long end = clock();

			cout << "n = " << n << "\n";
			cout << "size = " << size_3_17(n) << "\n";
			cout << "eps = " << eps << "\n";
			cout << "err = " << err << "\n";
			cout << "||z|| = " << z << "\n";
			cout << "||r|| = " << r << "\n";
			cout << "zeta = " << zeta << "\n";
			cout << "ro = " << ro << "\n";
			cout << "its = " << its << "\n";
			cout << "durability = " << (double)(end - start) / CLOCKS_PER_SEC << " s\n";
			cout << "\n\n\n";

			delete_vector(Ax, size_3_17(n));
			delete_vector(f, size_3_17(n));
			delete_vector(exact_sol, size_3_17(n));
			delete_vector(residual, size_3_17(n));
			delete_vector(error, size_3_17(n));
			delete_vector(sol, size_3_17(n));
		}
	}

	cout << "finish 1!!!\n\n\n\n";
	cout << "start 2!\n\n";
	for (int n = 100; n < 101; n++) {
		for (double eps = 1e-5; eps >= 1e-15; eps /= 10) {
			int its = 0;

			long start = clock();

			double* exact_sol = gen_exact_phi(n);
			double* f = gen_exact_f(n);

			double* sol = min_error_method_3_17(n, f, err, its, eps);

			double* error = vect_sub(sol, exact_sol, size_3_17(n));

			double* Ax = create_vector(size_3_17(n));
			mult_3_14(n, sol, Ax);

			double* residual = vect_sub(Ax, f, size_3_17(n));

			double z, r, zeta, ro;
			z = vect_norm(error, size_3_17(n));
			r = vect_norm(residual, size_3_17(n));
			zeta = z / vect_norm(exact_sol, size_3_17(n));
			ro = r / vect_norm(f, size_3_17(n));

			long end = clock();

			cout << "n = " << n << "\n";
			cout << "size = " << size_3_17(n) << "\n";
			cout << "eps = " << eps << "\n";
			cout << "err = " << err << "\n";
			cout << "||z|| = " << z << "\n";
			cout << "||r|| = " << r << "\n";
			cout << "zeta = " << zeta << "\n";
			cout << "ro = " << ro << "\n";
			cout << "its = " << its << "\n";
			cout << "durability = " << (double)(end - start) / CLOCKS_PER_SEC << " s\n";
			cout << "\n\n\n";

			delete_vector(Ax, size_3_17(n));
			delete_vector(f, size_3_17(n));
			delete_vector(exact_sol, size_3_17(n));
			delete_vector(residual, size_3_17(n));
			delete_vector(error, size_3_17(n));
			delete_vector(sol, size_3_17(n));
		}
	}

	cout << "finish 2!!!";

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

	//test();

	/*int n = 6;
	double* phi = gen_exact_phi(n);
	
	print_vector(phi, (n + 1) * (n + 1));

	double* f = gen_exact_f(n);

	double* phi_ = min_error_method_3_17(n, f);

	cout << "\n\n\n";

	print_vector(phi_, (n + 1) * (n + 1));

	delete_vector(f, (n + 1) * (n + 1));
	delete_vector(phi, (n + 1) * (n + 1));
	delete_vector(phi_, (n + 1) * (n + 1));
	*/

	test_2();

	for (int i = 0; i < 3; i++) {
		wait();
	}



	return 0;
}