// Copyright (c) 2018 Denis Skripov aka nizhikebinesi
//
// метод отражений с выбором ведущего столбца
#pragma once

#include <stdio.h>
#include <math.h>
#include <cstdlib>
#include <iostream>

using namespace std;

namespace numeric_methods {

	void wait() {
		char ch;
		//cin >> ch;
		scanf_s("%c", &ch);
	}

	void print_vector(double* vect, size_t n) {
		//cout.precision(3);
		for (size_t i = 0; i < n; i++) {
			//cout << vect[i] << " ";
			printf("%.3f ", vect[i]);
		}
	}

	//size_t* create_transpositions_vector(size_t n) {
	//	size_t* tmp = new size_t[n];
	//	return tmp;
	//}

	//void delete_transpositions_vector(size_t* vect, size_t n) {
	//	delete[] vect;
	//}



	//template <typename T> int sgn(T val) {
	//	return (T(0) < val) - (val < T(0));
	//}


	double* create_vector(size_t n) {
		double* tmp = new double[n];
		return tmp;
	}

	void zero_filling_vector(double* vect, size_t n) {
		for (size_t i = 0; i < n; i++) {
			vect[i] = 0;
		}
	}

	void random_int_filling_vector(double* vect, size_t n, size_t seed = 42, size_t divide_by = 10, size_t bias = 0) {
		srand(seed);
		for (size_t i = 0; i < n; i++) {
			vect[i] = rand() % divide_by + bias;
		}
	}

	void delete_vector(double* vect, size_t n) {
		delete[] vect;
	}

	double** create_matrix(size_t n) {
		double** tmp = new double*[n];
		for (size_t i = 0; i < n; i++) {
			tmp[i] = new double[n];
		}
		return tmp;
	}

	void zero_filling_matrix(double** matr, size_t n) {
		for (size_t i = 0; i < n; i++) {
			for (size_t j = 0; j < n; j++) {
				matr[i][j] = 0;
			}
		}
	}

	void random_int_filling_matrix(double** matr, size_t n, size_t seed = 42, size_t divide_by = 10, size_t bias = 0) {
		srand(seed);
		for (size_t i = 0; i < n; i++) {
			for (size_t j = 0; j < n; j++) {
				matr[i][j] = rand() % divide_by + bias;
			}
		}
	}

	void delete_matrix(double** matr, size_t n) {
		for (size_t i = 0; i < n; i++) {
			delete[] matr[i];
		}
		delete[] matr;
	}

	void print_matrix(double** matr, size_t n) {
		//cout.precision(3);
		for (size_t i = 0; i < n; i++) {
			for (size_t j = 0; j < n; j++) {
				//cout << matr[i][j] << " ";
				printf("%.3f ", matr[i][j]);
			}
			//cout << "\n";
			printf("\n");
		}
	}

	bool equal_doubles(double a, double b, double eps = 1e-9) {
		return fabs(a - b) < eps;
	}

	bool less_doubles(double a, double b, double eps = 1e-9) {
		return a < b && !equal_doubles(a, b, eps);
	}

	double infinity_norm(double** matr, size_t n) {
		double tmp(0);
		for (size_t i = 0; i < n; i++) {
			tmp += matr[i][0];
		}
		for (size_t j = 1; j < n; j++) {
			double t(0);
			for (size_t i = 0; i < n; i++) {
				t += matr[i][j];
			}
			if (less_doubles(tmp, t)) {
				tmp = t;
			}
		}
		return tmp;
	}

	double vector_length(double* a, size_t j, size_t size) {
		double tmp = 0.;
		for (size_t i = j; i < size; i++) {
			tmp += a[i] * a[i];
		}
		return sqrt(tmp);
	}

	/*void inverse_triangle_matrix(double** u, size_t size) {
		for (int i = size - 1; i > -1; i--) {
			double* u_i_n = create_vector(size - i);

			for (int l = i; l < size; l++) {
				u_i_n[l - i] = u[i][l];
			}

			u[i][i] = 1 / u_i_n[0];

			for (int j = i + 1; j < size; j++) {
				double tmp(0);
				for (size_t k = i + 1; k <= j; k++) {
					tmp += u_i_n[k - i] * u[k][j];
				}
				u[i][j] = tmp;
				u[i][j] /= -u_i_n[0];
			}

			delete_vector(u_i_n, size - i);
		}
	}*/

	size_t hand_filling(double**& a) {
		size_t size(0);
		scanf_s("%Id", &size);

		printf("%Id\n\n", size);

		a = create_matrix(size);

		for (size_t i = 0; i < size; i++) {
			for (size_t j = 0; j < size; j++) {
				double tmp(0);
				scanf_s("%lf", &tmp);

				//printf("%lf\n\n", tmp);
				a[i][j] = tmp;
			}
		}

		return size;
	}

	size_t first_filling(double**& a) {
		size_t size = 4;

		a = create_matrix(size);

		zero_filling_matrix(a, size);

		a[0][0] = 3;
		a[0][1] = 1;
		a[0][2] = 2;
		a[0][3] = 1;
		a[1][1] = 8;
		a[1][2] = 3;
		a[1][3] = 2;
		a[2][2] = 6;
		a[2][3] = 4;
		a[3][3] = 4;

		return size;
	}

	void swap_columns(double** a, size_t size, size_t i, size_t j) {
		for (int k = 0; k < size; k++) {
			double tmp = a[k][i];
			a[k][i] = a[k][j];
			a[k][j] = tmp;
		}
	}

	void swap_lines(double** a, size_t size, size_t i, size_t j) {
		double* tmp;
		tmp = a[i];
		a[i] = a[j];
		a[j] = tmp;
	}

	void delete_int_vector(int* t, size_t size) {
		delete[] t;
	}

	double scalar_product(double* a, double* b, size_t j, size_t size) {
		double tmp = 0.;
		for (size_t i = j; i < size; i++) {
			tmp += a[i] * b[i];
		}
		return tmp;
	}

	int index_for_swap(double** A, size_t size, int j) {
		int k = j;
		double sum(0);

		for (int i = j; i < size; i++) {
			sum += A[i][k] * A[i][k];
		}

		for (int k1 = j + 1; k1 < size; k1++) {
			double tmp(0);
			for (int i = j; i < size; i++) {
				tmp += A[i][k1] * A[i][k1];
			}

			if (less_doubles(sum, tmp)) {
				k = k1;
				sum = tmp;
			}
		}

		return k;
	}

	double** copy_matrix(double** a, size_t size) {
		double** tmp = create_matrix(size);
		for (size_t i = 0; i < size; i++) {
			for (size_t j = 0; j < size; j++) {
				tmp[i][j] = a[i][j];
			}
		}
		return tmp;
	}

	// todo! возможно ошибка здесь!
	void inverse_triangle_matrix(double** u, size_t size) {
		double *tmp = create_vector(size);
		for (int i = size - 1; i > -1; i--) {

			for (int j = i; j < size; j++) {
				tmp[j] = u[i][j];
			}

			u[i][i] = 1. / tmp[i];

			for (int j = i + 1; j < size; j++) {
				double t = 0.;
				for (size_t k = i + 1; k < j + 1; k++) {
					t -= tmp[k] * u[k][j];
				}
				u[i][j] = t / tmp[i];
			}

		}
		delete_vector(tmp, size);
	}

	void normal_w(double* a, double* w, size_t j, size_t size) {
		//for (size_t i = 0; i < j; i++) {
		//	w[i] = 0;
		//}
		for (size_t i = j; i < size; i++) {
			w[i] = -a[i];
		}

		double length_a = vector_length(w, j, size);

		w[j] += ((a[j] > 0) ? -length_a : length_a);

		double length = sqrt(2 * (length_a + fabs(a[j])) * length_a);
		for (size_t i = j; i < size; i++) {
			w[i] /= length;
		}
	}

	// still raw method
	// danger! todo this method!
	void just_reflection(double** A, size_t size) {
		double* diag = create_vector(size);
		for (int i = 0; i < size; i++) {
			diag[i] = A[i][i];
		}

		// column - лишнее
		double* column = create_vector(size);
		double* w = create_vector(size);

		for (int j = 0; j < size - 1; j++) {
			// берем редуцированный столбец, чтобы построить w
			for (int i = 0; i < j; i++) {
				column[i] = 0.;
				//w[i] = 0;
			}
			for (int i = j; i < size; i++) {
				column[i] = A[i][j];
				//w[i] = 0;
			}
		    normal_w(column, w, j, size);

			for (int i = j; i < size; i++) {
				// выделяем редуцированный столбец
				//for (int k = 0; k < j; k++) {
				//	column[k] = 0.;
				//}
				for (int k = j; k < size; k++) {
					column[k] = A[k][i];
				}
				double sc_prod = scalar_product(column, w, j, size);
				
				for (int k = j; k < size; k++) {
					A[k][i] -= 2. * sc_prod * w[k];
				}
			}

			diag[j] = A[j][j];
			for (int i = j; i < size; i++) {
				A[i][j] = w[i];
			}
		}
		
		// ошибка должна быть здесь!
		{
			double* w_diag = create_vector(size);

			for (int i = 0; i < size - 1; i++) {
				w_diag[i] = A[i][i];
				A[i][i] = diag[i];
			}

			// обратная треугольная
			inverse_triangle_matrix(A, size);

			for (int i = 0; i < size - 1; i++) {
				diag[i] = A[i][i];
				A[i][i] = w_diag[i];
			}

			delete_vector(w_diag, size);
		}

		//A[size - 1][size - 1] = diag[size - 1];
		for (int j = size - 2; j > -1; j--) {
			// берем редуцированный столбец w
			//for (int i = 0; i < j; i++) {
			//	w[i] = 0;
			//}
			for (int i = j; i < size; i++) {
				w[i] = A[i][j];
				A[i][j] = 0.;
			}
			A[j][j] = diag[j];

			for (int i = j; i < size; i++) {
				// выделяем редуцированную строку
				//for (int k = 0; k < j; k++) {
				//	column[k] = 0;
				//}
				for (int k = j; k < size; k++) {
					column[k] = A[i][k];
				}
				double sc_prod = scalar_product(column, w, j, size);

				for (int k = j; k < size; k++) {
					A[i][k] -= 2 * sc_prod * w[k];
				}
			}
		}
		

		//delete_vector(diag_w, size);
		delete_vector(column, size);
		delete_vector(w, size);
		delete_vector(diag, size);
	}

	// method for coursework
	// todo it!
	void reflection_method(double** A, size_t size) {
		int* transposes = new int[size];

		double* diag = create_vector(size);
		for (int i = 0; i < size; i++) {
			diag[i] = A[i][i];
		}

		delete_vector(diag, size);
		delete_int_vector(transposes, size);
	}

	double** matrix_mult(double** A, double** B, size_t size)
	{
		double** C = create_matrix(size);

		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				C[i][j] = 0;
				for (int k = 0; k < size; k++) {
					C[i][j] += A[i][k] * B[k][j];
				}
			}
		}

		return C;
	}
}