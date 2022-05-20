#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <utility> // std::pair

#include "../src/utils/umath.hpp"
#include "../src/utils/utimer.cpp"
#include "../src/utils/utils.hpp"
#include "../src/jacobi.hpp"
#include "../src/utils/uarg.hpp"

using namespace std;

// grid
vector<int> matrix_sizes = {128, 256, 512, 1024, 2056, 4096, 10000, 20000};
//vector<int> matrix_sizes = {500, 1000, 5000, 10000, 15000, 20000, 30000, 50000};
vector<int> n_workers = {1, 2, 4, 8, 12, 16, 20, 24, 28, 32};

// mean
int n_exec = 10;
float mean_seq, mean_th, mean_ff;
// Parameters
ulong max_iteration = 10000;
float tol           = 10e-6;
int verbose         = 0;

unsigned int seed = 20;
std::mt19937 gen(seed);
std::uniform_real_distribution<> dis(-2.0, 2.0);


inline static void 
init_linear_system(matrix_t &A, vector_t &b, const int n){
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

    cout << std::fixed;
    fstream fout;

    fout << std::fixed;
    fout << std::setprecision(2); 
  
    // opens an existing csv file or creates a new file.
    fout.open("data_parallel_analysis.csv", ios::out | ios::app);
    fout << "size" << ", "
         << "seq_time"   << ","
         << "th_time_1"  << ","
         << "th_time_2"  << ","
         << "th_time_4"  << ","
         << "th_time_8"  << ","
         << "th_time_12" << ","
         << "th_time_16" << ","
         << "th_time_20" << ","
         << "th_time_24" << ","
         << "th_time_28" << ","
         << "th_time_32" << ","
         << "ff_time_1"  << ","
         << "ff_time_2"  << ","
         << "ff_time_4"  << ","
         << "ff_time_8"  << ","
         << "th_time_12" << ","
         << "th_time_16" << ","
         << "th_time_20" << ","
         << "th_time_24" << ","
         << "th_time_28" << ","
         << "th_time_32" << ","
         << "\n";

    for (auto &n: matrix_sizes) {
        cout << "TEST: A dim:" << n << endl;
        // Data
        matrix_t A;
        vector_t b;
        init_linear_system(A, b, n);
        cout << "A size: "<< A.size()<< endl;

        if (!check_ddm(A)){
            cout << "ERROR: matrix is not ddm" << endl;
            return EXIT_FAILURE;
        }          

        mean_seq = 0;
        for (int i = 0; i < n_exec; i++){
            jacobi_seq(A, b, max_iteration, tol, verbose);
            mean_seq += jacobi_comp_time;
            //jacobi_th(A, b, max_iteration, tol, verbose, nw);
            //jacobi_ff(A, b, max_iteration, tol, verbose, nw);
        }
        mean_seq /= n_exec;
        cout << "seq mean time: "<< mean_seq << endl;
                // Insert the data to file
        fout << n << ", "
             << mean_seq << ",";

        for(auto &nw : n_workers){
            mean_th = 0;
            for (int i = 0; i < n_exec; i++){
                jacobi_th(A, b, max_iteration, tol, nw, verbose);
                mean_th += jacobi_comp_time;
            }
            mean_th /= n_exec;
            cout << "th_" << nw<< " mean time: "<< mean_th << endl;
            fout << mean_th << ",";
        }

        for(auto &nw : n_workers){
            mean_ff = 0;
            for (int i = 0; i < n_exec; i++){
                jacobi_ff(A, b, max_iteration, tol, nw, verbose);
                mean_ff += jacobi_comp_time;
            }
            mean_ff /= n_exec;
            cout << "ff_" << nw<< " mean time: "<< mean_ff << endl;
            fout << mean_ff << ",";
        }
        fout << "\n";
        cout << endl;
    }

    return EXIT_SUCCESS;
}


