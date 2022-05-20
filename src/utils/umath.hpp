#ifndef UMATH_H
#define UMATH_H

#include <random>
#include <vector>
#include <iostream>
#include <iomanip>      // std::setw
#include "utils.hpp"

typedef std::vector<std::vector<float> > matrix_t;
typedef std::vector<float> vector_t;

template<typename T>
T urand_t(T min, T max) {
    return (T(std::rand()) / T(RAND_MAX) * (max - min)) + min;
}

template<typename T>
void rnd_vector(const ulong size, std::vector<T> &v, const T min, const T max) {
    for (int i = 0; i < size; ++i) {
        v.emplace_back(urand_t(min, max));
    }
}

template<typename T>
void zeros_vector(const ulong size, std::vector<T> &v) {
    for (int i = 0; i < size; ++i) {
        v.emplace_back(0);
    }
}

template<typename T>
void ddm(const ulong size, std::vector<std::vector<T> > &matrix,
                                       const T min, const T max) {
    for (ulong i = 0; i < size; ++i) {
        T sum = T(0);
        std::vector<T> row_vec;
        for (ulong j = 0; j < size; ++j){
            T value = urand_t(min, max);
            sum += std::abs(value);
            row_vec.emplace_back(value);
        }
        sum -= abs(row_vec[i]);
        row_vec[i] = std::abs(row_vec[i]) + sum;
        matrix.emplace_back(row_vec);
    }
}

template<typename T>
bool check_ddm(std::vector<std::vector<T> > &matrix){
	int check_count = 0;
	// For each row ..
	for (int i = 0; i < matrix.size(); i++) {
		float row_sum = 0;
		// Summing the other row elements .. 
		for (int j = 0; j < matrix.size(); j++) {
			if (j != i) row_sum += abs(matrix[i][j]);
		}

		if (abs(matrix[i][i]) >= row_sum) {
			check_count++;
		}
	}
	return check_count == matrix.size();
}


template<typename T>
void cout_vec(const std::vector<T> &vec) {
    for (auto &value: vec) {
        std::cout << value << " ";
    }
    std::cout << std::endl;
}

template<typename T>
void cout_mat(const std::vector<std::vector<T> > &matrix) {
   for (auto &vector : matrix)
    {
        for (auto &x : vector)
        {
            std::cout << std::setw(12) << x ;
        }    
        std::cout << std::endl;
    }
}

#endif //UMATH_H