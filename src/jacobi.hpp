#ifndef JACOBI_H
#define JACOBI_H

#include <thread>
#include <vector>
#include <ff/parallel_for.hpp>

#include "utils/umath.hpp"
#include "utils/utils.hpp"

extern long jacobi_comp_time;

vector_t 
jacobi_seq(const matrix_t &A, const vector_t &b,
                             const uint64_t iter_max, const float tol, const int verbose);

vector_t 
jacobi_th(const matrix_t &A, const vector_t &b,
                             const uint64_t iter_max, const float tol, const uint64_t nw, const int verbose);

vector_t 
jacobi_ff(const matrix_t &A, const vector_t &b,
                             const uint64_t iter_max, const float tol, const uint64_t nw, const int verbose);

#endif // JACOBI_H
