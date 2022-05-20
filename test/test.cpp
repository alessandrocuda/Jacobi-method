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
ulong n             = DF_N;
ulong mode          = SEQ;
ulong nw            = 1;
ulong max_iteration = 10000;
float tol           = 10e-6;
int verbose         = 0;

// Data
matrix_t A;
vector_t b;
vector_t x;

unsigned int seed = 42;
std::mt19937 gen(seed);
std::uniform_real_distribution<> dis(1.0, 4.0);

inline static void 
cout_setup(const ulong n, const ulong mode, const ulong nw)
{
    cout << "Jacobi Test:" << endl;
    cout << "- generate a Linear system with " << n << " vars;" << endl;
    cout << "- solve with " << m.at(mode) << " with number of workers: " << nw << "." << endl;
}

inline static void 
init_linear_system(){
    for (ulong i = 0; i < n; ++i) {
        float sum = 0;
        vector_t row_vec;
        for (ulong j = 0; j < n; ++j){
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
    init_argv(argc, argv, n, mode, nw, seed, verbose);
    if (mode == SEQ) nw = 1;
    if (verbose) cout_setup(n, mode, nw);
    nw = min(nw, n);

    long seq_time;

    // Setup data
    if (verbose) cout << "Setup data...";
    init_linear_system();
    
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


