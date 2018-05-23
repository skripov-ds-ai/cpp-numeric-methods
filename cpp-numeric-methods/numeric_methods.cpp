// Copyright (c) 2018 Denis Skripov aka nizhikebinesi
//
// метод отражений с выбором ведущего столбца
#pragma once

#define _USE_MATH_DEFINES
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
			tmp += fabs(matr[i][0]);
		}
		for (size_t i = 1; i < n; i++) {
			double t(0);
			for (size_t j = 0; j < n; j++) {
				t += fabs(matr[i][j]);
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

	double** matrix_sub(double** A, double** B, size_t size) {
		double** C = create_matrix(size);

		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				C[i][j] = A[i][j] - B[i][j];
			}
		}

		return C;
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
		for (int i = size - 1; i >= 0; i--) {

			for (int j = i; j < size; j++) {
				tmp[j] = u[i][j];
			}

			u[i][i] = 1. / u[i][i];

			for (int j = i + 1; j < size; j++) {
				double t = 0.;
				for (int k = i + 1; k <= j; k++) {
					t += tmp[k] * u[k][j];
				}
				u[i][j] = - t / tmp[i];
			}

		}
		delete_vector(tmp, size);
	}

	void normal_w(double* a, double* w, size_t j, size_t size) {
		for (size_t i = j; i < size; i++) {
			w[i] = -a[i];
		}

		double length_a = vector_length(w, j, size);

		w[j] += ((a[j] > 0) ? -length_a : length_a);

		double length = sqrt(2 * length_a * (length_a + fabs(a[j])));
		for (size_t i = j; i < size; i++) {
			w[i] /= length;
		}
	}

	void yet_another_reflection(double** A, size_t size) {
		double* diag = create_vector(size);
		double* w = create_vector(size);

		double length, length_2, diag_el, dot_prod, length_v;

		for (int j = 0; j < size - 1; j++) {
			// вычисляем w
			length_2 = 0.;
			for (int i = j; i < size; i++) {
				w[i] = -A[i][j];
				length_2 += w[i] * w[i];
			}
			length = sqrt(length_2);
			diag_el = (A[j][j] > 0. ? -length : length);
			w[j] += diag_el;

			length_v = sqrt(2 * (length_2 + length * fabs(A[j][j])) );
			for (int i = j; i < size; i++) {
				w[i] /= length_v;
			}

			// применяем вектор w к стоблцам матрицы
			for (int i = j; i < size; i++) {
				dot_prod = 0.;
				for (int k = j; k < size; k++) {
					dot_prod += w[k] * A[k][i];
				}
				dot_prod += dot_prod;

				for (int k = j; k < size; k++) {
					A[k][i] -= dot_prod * w[k];
				}
			}
			
			// сохраняем w
			diag[j] = w[j];
			for (int i = j + 1; i < size; i++) {
				A[i][j] = w[i];
			}
		}

		// проверить диаг на ноль
		// обратная верхнетреугольная
		inverse_triangle_matrix(A, size);

		// обратное отражение
		for (int j = size - 2; j >= 0; j--) {
			w[j] = diag[j];
			for (int i = j + 1; i < size; i++) {
				w[i] = A[i][j];
				A[i][j] = 0.;
			}

			for (int i = 0; i < size; i++) {
				dot_prod = 0.;
				for (int k = j; k < size; k++) {
					dot_prod += w[k] * A[i][k];
				}
				dot_prod += dot_prod;

				for (int k = j; k < size; k++) {
					A[i][k] -= dot_prod * w[k];
				}
			}

		}

		delete_vector(w, size);
		delete_vector(diag, size);
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

	double** create_unit_matrix(size_t size) {
		double** E = create_matrix(size);

		zero_filling_matrix(E, size);

		for (size_t i = 0; i < size; i++) {
			E[i][i] = 1;
		}

		return E;
	}

	double* zero_phi_3_17(int n) {
		double* phi = create_vector((n + 1) * (n + 1));
		zero_filling_vector(phi, (n + 1) * (n + 1));
		
		return phi;
	}

	double* gen_exact_phi(int n) {
		double* phi = create_vector((n + 1) * (n + 1));
		double pi_2 = M_PI_2;

		for (int i = 0; i < n + 1; i++) {
			for (int j = 0; j < n + 1; j++) {
				phi[(n + 1) * i + j] = 0.;
			}
		}

		for (int i = 1; i < n; i++) {
			for (int j = 1; j < n; j++) {
				phi[(n + 1) * i + j] = sin(pi_2 * i) * sqrt(j);
				//cout << "\n" << sin(pi_2 * i) * sqrt(j) << "\n";
			}
		}

		return phi;
	}

	double* gen_exact_f(int n) {
		double* phi = gen_exact_phi(n);

		double* f = create_vector((n + 1) * (n + 1));
		for (int i = 0; i < n + 1; i++) {
			for (int j = 0; j < n + 1; j++) {
				f[(n + 1) * i + j] = 0.;
			}
		}

		for (int i = 1; i < n; i++) {
			for (int j = 1; j < n; j++) {
				f[(n + 1) * i + j] = -(
					phi[(n + 1) * i + (j + 1)] +
					phi[(n + 1) * i + (j - 1)] +
					phi[(n + 1) * (i + 1) + j] +
					phi[(n + 1) * (i - 1) + j] +
					-4 * phi[(n + 1) * i + j]
					) * ((double) n * n);
			}
		}


		delete_vector(phi, (n + 1) * (n + 1));

		return f;
	}

	void mult_3_14(int n, double* phi, double* f) {
		for (int i = 0; i < n + 1; i++) {
			for (int j = 0; j < n + 1; j++) {
				f[(n + 1) * i + j] = 0.;
			}
		}

		for (int i = 1; i < n; i++) {
			for (int j = 1; j < n; j++) {
				f[(n + 1) * i + j] = -(
					phi[(n + 1) * i + (j + 1)] +
					phi[(n + 1) * i + (j - 1)] +
					phi[(n + 1) * (i + 1) + j] +
					phi[(n + 1) * (i - 1) + j] +
					-4 * phi[(n + 1) * i + j]
					) * ((double) n * n);
			}
		}
	}

	// todo!
	double* min_error_method_3_17(int n, double* f, double eps=1e-9) {
		int size = (n + 1) * (n + 1);

		int it = 0;
		int maxit = 1000;
		
		double* tmp = create_vector(size);
		//double* tmp1 = create_vector(size);
		double* f_tmp = create_vector(size);
		//double* err = create_vector(size);
		double* delta = create_vector(size);
		double* tmp_phi = zero_phi_3_17(n);
		double* phi = zero_phi_3_17(n);

		double norm = 0.;
		double tau;
		
		do {
			mult_3_14(n, tmp_phi, f_tmp);

			// считаем невязку
			for (int i = 0; i < n + 1; i++) {
				for (int j = 0; j < n + 1; j++) {
					delta[(n + 1) * i + j] =
						f_tmp[(n + 1) * i + j] - f[(n + 1) * i + j];
				}
			}

			mult_3_14(n, delta, tmp);

			//cout << "\ndelta:\n";
			//print_vector(delta, size);

			tau = 0.;
			double denom = 0.;

			// считаем tau
			for (int i = 0; i < n + 1; i++) {
				for (int j = 0; j < n + 1; j++) {
					tau +=
						(
							delta[(n + 1) * i + j]
							*
							delta[(n + 1) * i + j]
						);
					denom += 
						(
							tmp[(n + 1) * i + j]
							*
							tmp[(n + 1) * i + j]
						);
				}
			}

			if (denom == 0) {
				it = -1;
				break;
			}

			//cout << "\ntau = " << tau << "\n";
			//cout << "\ndenom = " << denom << "\n";

			tau = tau / denom;

			//cout << "\ntmp:\n";
			//print_vector(tmp, size);
			//cout << "\nf:\n";
			//print_vector(f, size);

			for (int i = 0; i < n + 1; i++) {
				for (int j = 0; j < n + 1; j++) {
					phi[(n + 1) * i + j] =
						tmp_phi[(n + 1) * i + j] - tau * tmp[(n + 1) * i + j];
				}
			}

			norm = fabs(phi[0] - tmp_phi[0]);
			for (int i = 0; i < n + 1; i++) {
				for (int j = 0; j < n + 1; j++) {
					if (
						fabs(
							phi[(n + 1) * i + j] - tmp_phi[(n + 1) * i + j]
							) > norm
						) {
						norm = fabs(
							phi[(n + 1) * i + j] 
							-
							tmp_phi[(n + 1) * i + j]
						);
					}
					tmp_phi[(n + 1) * i + j] = phi[(n + 1) * i + j];
				}
			}

			//cout << "\nphi:\n";
			//print_vector(phi, size);
			//cout << "\nphi_tmp:\n";
			//print_vector(tmp_phi, size);


			it++;
		} while ((norm >= eps)/* & (it < maxit)*/);


		delete_vector(tmp, size);
		//delete_vector(tmp1, size);
		delete_vector(f_tmp, size);
		//delete_vector(err, size);
		delete_vector(delta, size);
		delete_vector(tmp_phi, size);
		//delete_vector(phi, size);
		return phi;
	}

	void matr_vect_mult(int size, double** A, double* v, double* res) {
		for (int i = 0; i < size; i++) {
			res[i] = 0.;
			for (int j = 0; j < size; j++) {
				res[i] += A[i][j] * v[j];
			}
		}
	}

	double* min_error_method(int n, double** A, double* f, double eps = 1e-9) {
		int it = 0;
		int maxit = 1000;

		double* tmp = create_vector(n);
		double* f_tmp = create_vector(n);
		double* delta = create_vector(n);
		double* tmp_phi = create_vector(n);
		double* phi = create_vector(n);

		double norm = 0.;
		double tau;

		do {
			matr_vect_mult(n, A, tmp_phi, f_tmp);

			// считаем невязку
			for (int i = 0; i < n; i++) {
				delta[i] =
						f_tmp[i] - f[i];
			}

			matr_vect_mult(n, A, delta, tmp);

			tau = 0.;
			double denom = 0.;

			// считаем tau
			for (int i = 0; i < n; i++) {
				tau +=
						(
							delta[i]
							*
							delta[i]
							);
				denom +=
						(
							tmp[i]
							*
							tmp[i]
							);
			}

			if (denom == 0) {
				it = -1;
				break;
			}

			tau = tau / denom;

			for (int i = 0; i < n; i++) {
				phi[i] =
					tmp_phi[i] - tau * tmp[i];
			}

			norm = fabs(phi[0] - tmp_phi[0]);
			for (int i = 0; i < n; i++) {
				if (
					fabs(
							phi[i] - tmp_phi[i]
						) > norm
					) {
						norm = fabs(
							phi[i]
							-
							tmp_phi[i]
						);
					}
					tmp_phi[i] = phi[i];		
			}

			it++;
		} while ((norm >= eps)/* & (it < maxit)*/);


		delete_vector(tmp, n);
		delete_vector(f_tmp, n);
		delete_vector(delta, n);
		delete_vector(tmp_phi, n);
		return phi;
	}
}