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

    inline void init_solutions(const ulong n){
        x.clear();
        x_next.clear();
        zeros_vector(n, x);
        zeros_vector(n, x_next);
    }

    inline void init_stop_condition(){
        error = FLT_MAX;
        //stopping_count = 0;
    }

    inline void compute_error(const ulong n, const float tol){
        float local_error = 0;
        // has been vectorized, not present with -fopt-info-vector-missed
        for (ulong i = 0; i < n; ++i){
            local_error += std::abs(x[i] - x_next[i]);
        }
        error = local_error/n;
    }

    inline void compute_x_next(const matrix_t &A, const vector_t &b, const int i){
        float sigma = 0;
        // has been vectorized, not present with -fopt-info-vector-missed
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
    long norm_swap_time = 0;
    {
        utimer timer("Jacobi main loop", &jacobi_comp_time, verbose);
        for (ulong iter = 0; iter < iter_max && error > tol; ++iter) {
            //compute x
            for (int i = 0; i < n; ++i) {
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
partial_jacobi(const ulong th_id, const ulong n_chunks, const ulong nw,
               const matrix_t &A, const vector_t &b, barrier &spin_barrier,
               const ulong iter_max, const float tol, const int verbose){

    //compute ranges
    ulong n = x.size();
    ulong start = th_id * n_chunks;
    ulong end = (th_id != nw - 1 ? start + n_chunks : n) - 1;
    
    if (verbose > 2){
        const std::lock_guard<std::mutex> lock(cout_mutex);
        std::cout << "TH_"<< th_id << \
                    ": (" << start <<", "<< end<< ")" << \
                    "#rows: " << (end - start + 1)<< std::endl;
    }


    //Start partial Jacobi method
    for (ulong iter = 0; iter < iter_max && error > tol; ++iter) {
        //calculate partizal x
        for (int i = start; i <= end; ++i) {
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
          const ulong iter_max, const float tol, const ulong nw, const int verbose){

    //initialize solution with zeros
    ulong n = A.size();
    init_solutions(n); 
    init_stop_condition();

    //initialize barrier and vector of threads 
    barrier spin_barrier(nw);
    std::vector<std::thread> threads;
    //compute #chunks 
    ulong n_chunks = n / nw;

    {
        utimer timer("Jacobi main loop", &jacobi_comp_time, verbose);
        for (ulong i = 0; i < nw; i++) {
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
          const ulong iter_max, const float tol, const ulong nw, const int verbose){
    
    //initialize solution with zeros
    ulong n = A.size();
    init_solutions(n); 
    init_stop_condition();

    {
        utimer timer("Jacobi main loop", &jacobi_comp_time, verbose);
        ff::ParallelFor pf(nw, true, true);
        // In case of static scheduling (chunk <= 0), the scheduler thread is never started.
        // pf.disableScheduler();

        for (ulong iter = 0; iter < iter_max && error > tol; ++iter) {
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

// partial_jacobi: MAP+reduce (partial matrix-vector moltiplication + local error computation)
/*
void 
partial_jacobi(const ulong th_id, const ulong n_chunks, const ulong nw,
               const matrix_t &A, const vector_t &b, barrier &spin_barrier,
               const ulong iter_max, const float tol, const int verbose){

    //compute ranges
    ulong n = x.size();
    ulong start = th_id * n_chunks;
    ulong end = (th_id != nw - 1 ? start + n_chunks : n) - 1;
    float local_error;    
    if (verbose > 1){
        const std::lock_guard<std::mutex> lock(cout_mutex);
        std::cout << "TH_"<< th_id << \
                    ": (" << start <<", "<< end<< ")" << \
                    "#rows: " << (end - start + 1)<< std::endl;
    }


    //Start partial Jacobi method
    for (ulong iter = 0; iter < iter_max && error > tol; ++iter) {
        //calculate partizal x
	local_error = 0;
        for (int i = start; i <= end; ++i) {
            compute_x_next(A, b, i);
	        local_error += std::abs(x[i] - x_next[i]);
        }

	    g_r_error[th_id] = local_error;
        
        spin_barrier.busywait();

        //#pragma omp once: error computed only by one threads 
        if (!flag.test_and_set()) {
          //compute_error(n, tol);
	        error = 0;
	        for (int i = 0; i < g_r_error.size(); ++i)
		        error += g_r_error[i];
	        error/=n;
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
*/

// FastFlow implementation with MAP+reduce
/* 
vector_t 
jacobi_ff(const matrix_t &A, const vector_t &b,
          const ulong iter_max, const float tol, const ulong nw, const int verbose){
    
    //initialize solution with zeros
    ulong n = A.size();
    init_solutions(n); 
    init_stop_condition();

    {
        utimer timer("Jacobi main loop", &jacobi_comp_time, verbose);
        ff::ParallelForReduce<double> pf(nw, true, true);
        // In case of static scheduling (chunk <= 0), the scheduler thread is never started.
        // pf.disableScheduler();
	    double lerror;
        for (ulong iter = 0; iter < iter_max && error > tol; ++iter) {
            lerror = 0.0;
		    pf.parallel_reduce(lerror, 0.0, 
                                0, n, 
                                1, 0,
                                [&](const long i, double& lerror) {
                                    compute_x_next(A, b, i);
				                    lerror += std::abs(x[i] - x_next[i]); 
                                },
			                    [](double& s, const double& d) { s+=d;},
                                nw);
	        error = lerror/n;
            //compute_error(n, tol);
            if (verbose > 1) std::cout << "Iter: "<< iter << std::setw(10) << "Error: " << error << std::endl;
            swap(x, x_next);    
        }
    }
    return x;

}
*/

/*
// At high value of n the last computation may slow down the whole outer jacobi iteration.

vector_t 
jacobi_ff_v2(const matrix_t &A, const vector_t &b,
          const ulong iter_max, const float tol, const ulong nw, const int verbose){
    
    //initialize solution with zeros
    ulong n = A.size();
    init_solutions(n); 
    init_stop_condition();

    ff::ParallelFor pf(nw-1, true, true);

    {
        utimer timer("Jacobi main loop", &jacobi_comp_time, verbose);
        int end_internal_range = n/nw;
        for (ulong iter = 0; iter < iter_max && error > tol; ++iter) {
            pf.parallel_for(end_internal_range, n, 1, 0,
                            [&](const long i) {
                                compute_x_next(A, b, i);
                            });
            for (size_t i = 0; i < end_internal_range; i++){
                compute_x_next(A, b, i);
            }
            compute_error(n, tol);
            if (verbose > 1) std::cout << "Iter: "<< iter << std::setw(10) << "Error: " << error << std::endl;
            std::swap(x, x_next);    
        }
    }
    return x;

}
*/