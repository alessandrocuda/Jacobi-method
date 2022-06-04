```bash
# Jemalloc
cd git clone https://github.com/jemalloc/jemalloc.git
cd jmalloc
./autogen.sh
make
cd ../test

# use JMALLOC to reduce the numbero of lock/unlock to the heap
LD_PRELOAD=../jemalloc/lib/libjemalloc.so.2 ./bin/test -n 10000 -m th -w 16 -t 10e-6 -s 20 -v 1
```

```c++
// partial_jacobi: 
//    - MAP+reduce (partial matrix-vector moltiplication  local error computation) 
//    - std::atomic_flag flag to execute the global error only once.
//    - require double barrier (less efficient)

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

	    th_error[th_id] = local_error;
        
        spin_barrier.busywait();

        //#pragma omp once: error computed only by one threads 
        if (!flag.test_and_set()) {
          //compute_error(n, tol);
	        error = 0;
	        for (int i = 0; i < th_error.size(); ++i)
		        error += th_error[i];
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

```

```c++
// FastFlow implementation with MAP+reduce
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
                                    lerror /=n;
                                },
			                    [](double& s, const double& d) { s+=d;},
                                nw);
            //compute_error(n, tol);
            if (verbose > 1) std::cout << "Iter: "<< iter << std::setw(10) << "Error: " << error << std::endl;
            swap(x, x_next);    
        }
    }
    return x;

}

```c++
// fastflow: leave last chuck to be computed by the main thread, to remove a thread. 
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
```