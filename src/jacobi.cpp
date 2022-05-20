#include "jacobi.hpp"
#include "barrier/barrier.hpp"
#include <mutex>
#include "utils/utimer.cpp"
#include <typeinfo>

namespace 
{  
    std::mutex cout_mutex;  // protects cout
    std::atomic_flag flag;  // compute error only once in jacobi_th method

    int stopping_count;
    float error;

    vector_t x_next;
    vector_t x;

    inline void init_solutions(const ulong n){
        x.clear();
        x_next.clear();
        zeros_vector(n, x);
        zeros_vector(n, x_next);
    }

    inline void init_stop_condition(){
        error = 0;
        stopping_count = 0;
    }

    inline void compute_stop_condition(const ulong n, const float tol){
        error = 0;
        stopping_count = 0;
        float absolute_error;
        for (ulong i = 0; i < n; ++i){
            absolute_error = std::abs(x[i] - x_next[i]);
            if (absolute_error <= tol) stopping_count++;
            error += absolute_error;
        }
        error /= n;
    }

    inline void compute_x_next(const matrix_t &A, const vector_t &b, const int i){
        float sigma = 0;
        for (int j = 0; j < A.size(); ++j) {
            sigma += (x[j] * A[i][j]);
        }
        x_next[i] = (-sigma + b[i] + x[i]*A[i][i]) / A[i][i];
    }
}

long jacobi_comp_time;

vector_t 
jacobi_seq(const matrix_t &A, const vector_t &b,
           const ulong iter_max, const float tol, const int verbose) {

    ulong n = A.size();
    //initialize solution with zeros
    init_solutions(n);
    init_stop_condition();

    //Start Jacobi method
    //int norm_swap_time = 0;
    {
        utimer timer("Jacobi main loop", &jacobi_comp_time, verbose);
        for (ulong iter = 0; iter < iter_max && stopping_count != n; ++iter) {
            //calculate x
            for (int i = 0; i < n; ++i) {
                compute_x_next(A, b, i);
            }
            //compute the error
            //START(start_norm_swap_time);

            compute_stop_condition(n, tol);
            if (verbose > 1) std::cout << "Iter: "<< iter << std::setw(10) << "Error: " << error << std::endl;
            swap(x, x_next);

            // STOP(start_norm_swap_time, elapsed_norm_swap_time);  
            //norm_swap_time += elapsed_norm_swap_time;
        }
    }
    return x;
}


void 
partial_jacobi(const ulong th_id, const ulong start, const ulong end,
               const matrix_t &A, const vector_t &b, barrier &spin_barrier,
               const ulong iter_max, const float tol, const int verbose){

    if (verbose > 1){
        const std::lock_guard<std::mutex> lock(cout_mutex);
        std::cout << "TH_"<< th_id << \
                    ": (" << start <<", "<< end<< ")" << \
                    "#rows: " << (end - start + 1)<< std::endl;
    }
    
    int n = x.size();
    init_stop_condition();


    //Start partial Jacobi method
    for (ulong iter = 0; iter < iter_max && stopping_count != n; ++iter) {
        //calculate partizal x
        for (int i = start; i <= end; ++i) {
            compute_x_next(A, b, i);
        }
        
        spin_barrier.busywait();

        //#pragma omp once: error computed only by one threads 
        if (!flag.test_and_set()) {
            compute_stop_condition(n, tol);
            if (verbose > 1) std::cout << "TH:"<< th_id<< " - iter: "<< iter << " - Error: " << error << std::endl;
            swap(x, x_next);    
        }

        spin_barrier.busywait();
        flag.clear();
    }
    
    if (verbose > 1){
        const std::lock_guard<std::mutex> lock(cout_mutex);
        std::cout << "TH_"<< th_id << ": exit"<< std::endl;
    }
    
}

vector_t 
jacobi_th(const matrix_t &A, const vector_t &b,
          const ulong iter_max, const float tol, const ulong nw, const int verbose){

    //initialize solution with zeros
    ulong n = A.size();
    init_solutions(n); 
    //initialize barrier and vector of threads 
    barrier spin_barrier(nw);
    std::vector<std::thread> threads;
    //compute #chunks 
    ulong n_chunks = n / nw;

    {
        utimer timer("Jacobi main loop", &jacobi_comp_time, verbose);
        for (ulong i = 0; i < nw; i++) {
            //compute ranges and pass them to the thread
            ulong start = i * n_chunks;
            ulong end = (i != nw - 1 ? start + n_chunks : n) - 1;
            threads.emplace_back(std::thread(partial_jacobi, i, start, end,
                                                std::ref(A), std::ref(b), 
                                                std::ref(spin_barrier),
                                                iter_max, tol, verbose));
        }

        for (auto &th: threads)
            th.join();    
    }


    return x;
}


vector_t 
jacobi_ff(const matrix_t &A, const vector_t &b,
          const ulong iter_max, const float tol, const ulong nw, const int verbose){
    
    //initialize solution with zeros
    ulong n = A.size();
    init_solutions(n); 
    init_stop_condition();

    ff::ParallelFor pf(nw, true, true);

    {
        utimer timer("Jacobi main loop", &jacobi_comp_time, verbose);
        for (ulong iter = 0; iter < iter_max && stopping_count != n; ++iter) {
            pf.parallel_for(0, n, 1, 0,
                            [&](const long i) {
                                compute_x_next(A, b, i);
                            },
                            nw);
            compute_stop_condition(n, tol);
            if (verbose > 1) std::cout << "Iter: "<< iter << std::setw(10) << "Error: " << error << std::endl;
            swap(x, x_next);    
        }
    }
    return x;

}


