#include <string>
#include <iostream>

#include "../src/utils/umath.hpp"
#include "../src/utils/utimer.cpp"
#include "../src/utils/uarg.hpp"
#include "../src/utils/utils.hpp"
#include "../src/jacobi.hpp"
#include "../src/utils/uarg.hpp"

using namespace std;

// Parameters
uint64_t n             = DF_N;
uint64_t mode          = SEQ;
uint64_t nw            = 1;
uint64_t max_iteration = 10000;
float tol           = 10e-7;
unsigned int seed   = 20;
std::mt19937 gen(seed);
float l_range       = -1.0;
float r_range       =  8.0;
int verbose         = 0;

// Data
matrix_t A;
vector_t b;
vector_t x;



inline static void 
cout_setup(const uint64_t n, const uint64_t mode, const uint64_t nw)
{
    cout << "Jacobi Test:" << endl;
    cout << "- generate a Linear system with " << n << " vars;" << endl;
    cout << "- solve with " << m.at(mode) << " with number of workers: " << nw << "." << endl;
}

inline static void 
init_linear_system(){
    // not present with -fopt-info-vector-missed 
    std::uniform_real_distribution<> dis(l_range, r_range);//dis(-1.0, 8.0);
    for (size_t i = 0; i < n; ++i) {
        float sum = 0;
        vector_t row_vec;
        for (size_t j = 0; j < n; ++j){
            float value = dis(gen);
            sum += abs(value);
            row_vec.emplace_back(value);
        }
        sum -= abs(row_vec[i]);
        row_vec[i] = abs(row_vec[i]) + sum;
        A.emplace_back(row_vec);
        b.emplace_back(dis(gen));
    }
}

int
main(int argc, char *const argv[]){

    // Setup parameters
    init_argv(argc, argv, n, mode, nw, seed, l_range, r_range, tol, verbose);
    if (mode == SEQ) nw = 1;
    if (verbose) cout_setup(n, mode, nw);
    nw = min(nw, n);

    long seq_time;

    // Setup data
    if (verbose) cout << "Setup data...";

    {
        utimer timer("init linear system");
        init_linear_system();
    }
    
    //static const float range = 4.f;
    // ddm(n, A, range, range);
    // rnd_vector(n, b, range, range);
    if (!check_ddm(A)){
        cout << "ERROR: matrix is not ddm" << endl;
        return EXIT_FAILURE;
    } 

    if(verbose and n < 20){
        cout << endl;
        cout_mat(A);
        cout << endl;
        cout_vec(b);
    }

    if (verbose) cout << "solve Ax=b..."<< endl;
    //Select method
    switch (mode)
    {
    case SEQ:
        {
        utimer timer("jacobi_seq", &seq_time);
        x = jacobi_seq(A, b, max_iteration, tol, verbose);
        }
        break;
    case TH:
        {
        utimer timer("jacobi_th", &seq_time);
        x = jacobi_th(A, b, max_iteration, tol, nw, verbose);
        }
        break;
    case FF:
        {
        utimer timer("jacobi_ff", &seq_time);
        x = jacobi_ff(A, b, max_iteration, tol, nw, verbose);
        }
        break;
    default:
        break;
    }

    if(verbose and n < 20) cout_vec(x);
    return EXIT_SUCCESS;
}


