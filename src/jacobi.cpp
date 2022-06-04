#include "jacobi.hpp"
#include "barrier/barrier.hpp"
#include <mutex>
#include "utils/utimer.cpp"
#include <typeinfo>
#include <cfloat>


namespace 
{  
    std::mutex cout_mutex;  // protects cout
    std::atomic_flag flag;  // compute error only once in jacobi_th method

    int stopping_count;
    float error;

    vector_t x_next;
    vector_t x;

    inline void init_solutions(const uint64_t n){
        x.clear();
        x_next.clear();
        zeros_vector(n, x);
        zeros_vector(n, x_next);
    }

    inline void init_stop_condition(){
        error = FLT_MAX;
        //stopping_count = 0;
    }

    inline void compute_error(const uint64_t n, const float tol){
        float local_error = 0;
        // has been vectorized, not present with -fopt-info-vector-missed
        for (size_t i = 0; i < n; ++i){
            local_error += std::abs(x[i] - x_next[i]);
        }
        error = local_error/n;
    }

    inline void compute_x_next(const matrix_t &A, const vector_t &b, const int i){
        float sigma = 0;
        // has been vectorized, not present with -fopt-info-vector-missed
        for (size_t j = 0; j < A.size(); ++j) {
            sigma += (x[j] * A[i][j]);
        }
        x_next[i] = (-sigma + b[i] + x[i]*A[i][i]) / A[i][i];
    }
}

long jacobi_comp_time;

vector_t 
jacobi_seq(const matrix_t &A, const vector_t &b,
           const uint64_t iter_max, const float tol, const int verbose) {

    uint64_t n = A.size();
    //initialize solution with zeros
    init_solutions(n);
    init_stop_condition();

    //Start Jacobi method
    long norm_swap_time = 0;
    {
        utimer timer("Jacobi main loop", &jacobi_comp_time, verbose);
        for (size_t iter = 0; iter < iter_max && error > tol; ++iter) {
            //compute x
            for (size_t i = 0; i < n; ++i) {
                compute_x_next(A, b, i);
            }
            //compute the error
            START(start_norm_swap_time);
            compute_error(n, tol);
            if (verbose > 1) std::cout << "Iter: "<< iter << std::setw(10) << "Error: " << error << std::endl;
            std::swap(x, x_next);

            STOP(start_norm_swap_time, elapsed_norm_swap_time);  
            norm_swap_time += elapsed_norm_swap_time;
        }
    }
    if (verbose > 0)
    std::cout << "Jacobi time: " << jacobi_comp_time << " norm time: " << norm_swap_time << std::endl;
    return x;
}


void 
partial_jacobi(const uint64_t th_id, const uint64_t n_chunks, const uint64_t nw,
               const matrix_t &A, const vector_t &b, barrier &spin_barrier,
               const uint64_t iter_max, const float tol, const int verbose){

    //compute ranges
    uint64_t n = x.size();
    uint64_t start = th_id * n_chunks;
    uint64_t end = (th_id != nw - 1 ? start + n_chunks : n) - 1;
    
    if (verbose > 2){
        const std::lock_guard<std::mutex> lock(cout_mutex);
        std::cout << "TH_"<< th_id << \
                    ": (" << start <<", "<< end<< ")" << \
                    "#rows: " << (end - start + 1)<< std::endl;
    }


    //Start partial Jacobi method
    for (size_t iter = 0; iter < iter_max && error > tol; ++iter) {
        //calculate partizal x
        for (size_t i = start; i <= end; ++i) {
            compute_x_next(A, b, i);
        }
        //barrier + #pragma omp once: error computed only by one threads 
        spin_barrier.busywait([&] {
            compute_error(n, tol);
            if (verbose > 1) std::cout << "TH:"<< th_id<< " - iter: "<< iter << " - Error: " << error << std::endl;
            std::swap(x, x_next);    
        });
    }
    
    if (verbose > 2){
        const std::lock_guard<std::mutex> lock(cout_mutex);
        std::cout << "TH_"<< th_id << ": exit"<< std::endl;
    }
    
}

vector_t 
jacobi_th(const matrix_t &A, const vector_t &b,
          const uint64_t iter_max, const float tol, const uint64_t nw, const int verbose){

    //initialize solution with zeros
    uint64_t n = A.size();
    init_solutions(n); 
    init_stop_condition();

    //initialize barrier and vector of threads 
    barrier spin_barrier(nw);
    std::vector<std::thread> threads;
    //compute #chunks 
    uint64_t n_chunks = n / nw;

    {
        utimer timer("Jacobi main loop", &jacobi_comp_time, verbose);
        for (size_t i = 0; i < nw; i++) {
            threads.emplace_back(std::thread(partial_jacobi, i, n_chunks, nw,
                                                std::ref(A), std::ref(b), 
                                                std::ref(spin_barrier),
                                                iter_max, tol, verbose));
        }

        for (auto &th: threads)
            th.join();    
    }

    if (verbose > 0)
    std::cout << "Jacobi time: " << jacobi_comp_time << std::endl;
    
    return x;
}


vector_t 
jacobi_ff(const matrix_t &A, const vector_t &b,
          const uint64_t iter_max, const float tol, const uint64_t nw, const int verbose){
    
    //initialize solution with zeros
    uint64_t n = A.size();
    init_solutions(n); 
    init_stop_condition();

    {
        utimer timer("Jacobi main loop", &jacobi_comp_time, verbose);
        ff::ParallelFor pf(nw, true, true);
        // In case of static scheduling (chunk <= 0), the scheduler thread is never started.
        // pf.disableScheduler();

        for (size_t iter = 0; iter < iter_max && error > tol; ++iter) {
            pf.parallel_for(0, n, 1, 0,
                            [&](const long i) {
                                compute_x_next(A, b, i);
                            },
                            nw);
            compute_error(n, tol);
            if (verbose > 1) std::cout << "Iter: "<< iter << std::setw(10) << "Error: " << error << std::endl;
            std::swap(x, x_next);    
        }
    }
    if (verbose > 0)
    std::cout << "Jacobi time: " << jacobi_comp_time << std::endl;
    
    return x;

}