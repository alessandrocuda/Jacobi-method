# 


# Jacobi-method
[![alessandrocuda - Jacobi-method](https://img.shields.io/static/v1?label=alessandrocuda&message=Jacobi-method&color=blue&logo=github)](https://github.com/alessandrocuda/Jacobi-method "Go to GitHub repo")
[![stars - Jacobi-method](https://img.shields.io/github/stars/alessandrocuda/Jacobi-method?style=social)](https://github.com/alessandrocuda/Jacobi-method)
[![forks - Jacobi-method](https://img.shields.io/github/forks/alessandrocuda/Jacobi-method?style=social)](https://github.com/alessandrocuda/Jacobi-method)


the Jacobi method is an iterative algorithm for determining the solutions of a strictly diagonally dominant system of linear equations.

## Details
TBW

## Table of Contents 
- [Build](#usage)
- [Contributing](#contributing)
- [Contact](#contact)
- [License](#license)

## Build
```bash
# clone FastFlow Lib in a folder
git clone https://github.com/fastflow/fastflow.git
# Then set the FF env variable
export FF_ROOT=$PWD/fastflow/
# NOTE: remember to run the Bash script mapping_string.sh under fastflow/ff/
cd fastflow/ff && ./mapping_string.sh && cd ../..
# clone the jacobi repo in a folder
git clone https://github.com/alessandrocuda/Jacobi-method.git
cd Jacobi-method
# compile all with the makefile
make -j
# or you can choose to compile the desired executable file
make bin/test
make bin/parallel_analysis
```

After the compilation two bin files are generated: `bin/test` and `bin/parallel_analysis`. Here, all the command line parameters are reported:

1. `./bin/test [-n N] [-m mode] -w [nw] [-s seed] [-l lv] [-r rv] [-t tol] [-h] [-v]`
    The -h flag will print the helper:
    ```
    Parameters related on problem dimension, method and parallel degree:
        -n   N           matrix dimension NxN, N >1, default N=10
        -m   mode        implementation mode = [seq, th, ff], default mode=seq
                        seq: sequential, th: threads, ff: fastflow
        -w   nw          number of workers/threads (nw) used in the 
                        paralell computation, default nw=1
        -t   tol         the stop condition error tolerance, default tol=10e-7
        
    Parameters used to generate the Linear system: set seed, draws from a uniform 
    distribution with ranges U[lv, rv], then force it strictly diagonal dominant:
        -s   seed        default seed= 20
        -l   lv          default lv  = -1.0
        -r   rv          default rv  = -8.0            
        -h               print this help and exit
        -v               verbose level, 0: None, 1:print problem, results and times, 
                        2: print 1 and every iteration of the algorithm  
    ```

2. `./bin/parallel_analysis file_name}`:

    It performs a statistical analysis by varying the magnitude of the problem, the method and the number of threads. For each iteration, the mean is computed and saved in a user-defined file `file_name`. Ranges are defined as follow:

        matrix_sizes = {500, 1000, 5000, 10000, 15000, 20000, 30000, 50000};
        n_workers    = {1, 2, 4, 8, 12, 16, 20, 24, 28, 32};
    


## Contributing
 
1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request :D

<!-- CONTACT -->
## Contact

Alessandro Cudazzo - [@alessandrocuda](https://twitter.com/alessandrocuda) - alessandro@cudazzo.com


<!-- LICENSE -->
## License
[![License](http://img.shields.io/:license-mit-blue.svg?style=flat-square)](http://badges.mit-license.org)

This library is free software; you can redistribute it and/or modify it under
the terms of the MIT license.

- **[MIT license](LICENSE)**
- Copyright 2019 ©  <a href="https://alessandrocudazzo.it" target="_blank">Alessandro Cudazzo</a> 
