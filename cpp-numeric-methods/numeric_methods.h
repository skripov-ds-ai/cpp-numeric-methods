#pragma once

namespace numeric_methods {

	void wait();

	void print_vector(double* vect, size_t n);

	template <typename T> int sgn(T val) {
		return (T(0) < val) - (val < T(0));
	}

	double* create_vector(size_t n);

	void zero_filling_vector(double* vect, size_t n);

	void random_int_filling_vector(double* vect, size_t n, size_t seed = 42, size_t divide_by = 10, size_t bias = 0);

	void delete_vector(double* vect, size_t n);

	double** create_matrix(size_t n);

	void zero_filling_matrix(double** matr, size_t n);

	void random_int_filling_matrix(double** matr, size_t n, size_t seed = 42, size_t divide_by = 10, size_t bias = 0);

	void delete_matrix(double** matr, size_t n);

	void print_matrix(double** matr, size_t n);

	bool equal_doubles(double a, double b, double eps = 1e-9);

	bool less_doubles(double a, double b, double eps = 1e-9);

	double infinity_norm(double** matr, size_t n);

	double vector_length(double* a, size_t j, size_t size);

	void inverse_triangle_matrix(double** u, size_t size);


	void normal_w(double* a, double*w, size_t j, size_t size);

	size_t hand_filling(double**& a);

	size_t first_filling(double**& a);

	void swap_columns(double** a, size_t size, size_t i, size_t j);

	void swap_lines(double** a, size_t size, size_t i, size_t j);

	void delete_int_vector(int* t, size_t size);

	double scalar_product(double* a, double* b, size_t j, size_t size);

	int index_for_swap(double** A, size_t size, int j);

	double** copy_matrix(double** a, size_t size);

	void yet_another_reflection(double** A, size_t size);

	double** create_unit_matrix(size_t size);

	double** matrix_sub(double** A, double** B, size_t size);

	double** matrix_mult(double** A, double** B, size_t size);

	void mult_3_14(int n, double* phi, double* f);

	double* gen_exact_f(int n);

	double* gen_exact_phi(int n);

	double* zero_phi_3_17(int n);

	double* min_error_method_3_17(int n, double* f, double eps = 1e-9);
}