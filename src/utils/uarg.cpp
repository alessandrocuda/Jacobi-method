#include "uarg.hpp"

static void
cout_usage()
{
    fputs("Jacobi Test: [OPTIONS]\n", stdout);
	fputs("Parameters related on problem dimension, method and parallel degree:",stdout);
    fputs("-n   N           matrix dimension NxN, N >1, default N=10 \n", stdout);
    fputs("-m   mode        implementation mode = [seq, th, ff], default mode=seq \n", stdout);
    fputs("                 seq: sequential, th: threads, ff: fastflow \n", stdout);
    fputs("-w   nw          number of workers/threads (nw) used in the paralell computation, default nw=1\n", stdout);
    fputs("-t   tol         the stop condition error tolerance, default tol=10e-7 \n", stdout);
    fputs("\n", stdout);
	fputs("Parameters related on how the linear system is generated: set a certain seed, \n", stdout);
	fputs("draws from a uniform distribution with ranges A_ij in U[lv, rv], \n",stdout);
	fputs("and  made it strictly diagonal dominant: \n ",stdout);
    fputs("-s   seed        default seed= 20\n", stdout);
    fputs("-l   lv          default lv  = -1.0\n", stdout);
    fputs("-r   rv          default rv  = -8.0\n", stdout);
    fputs("-h               print this help and exit\n", stdout);
    fputs("-v               verbose level, 0: None, 1:print problem, results and times, \n", stdout);
}

void init_argv(const int argc, char *const argv[], 
					uint64_t &n, uint64_t &mode, uint64_t &nw, 
					unsigned int &seed, float &l_range, float &r_range, 
					float &tol, int &verbose) {
	int 	c; 					// temp arg for getopt, GNU c lib std

	while ((c = getopt(argc, argv, "n:m:w:s:l:r:t:hv:")) != -1){
		switch (c){
			case 'n':
				if (isdigit(*optarg)){
					n = (uint64_t) atoi(optarg);
				}
				break;
			case 'm':
				if (strcmp(optarg, SEQUENTIAL) == 0) mode = SEQ;
				else if (strcmp(optarg, THREADS) == 0) mode = TH;
				else if (strcmp(optarg, FASTFLOW) == 0) mode = FF;
				else{
					fputs("Invalid method, set to sequential.\n", stdout);
					mode = SEQ;
				}
				break;
			case 'w':
				if (isdigit(*optarg)){
					nw = (uint64_t) atoi(optarg);
				}
				break;
			case 's':
				if (isdigit(*optarg)){
					seed = (uint64_t) atoi(optarg);
				}
				break;
			case 'l':
				l_range =  atof(optarg);
				break;
			case 'r':
				r_range =  atof(optarg);
				break;
			case 't':
				tol =  atof(optarg);
				break;
            case 'v':
				if (isdigit(*optarg)){
					verbose = (uint64_t) atoi(optarg);
				}
                break;
			case 'h':
				cout_usage();
				exit(EXIT_SUCCESS);
			case '?':
				if (isprint (optopt)){
				  	fprintf ( stderr, 
				  			  "Unknown option `-%c'.\n", optopt);
					 cout_usage();
				}
				else{
					//char is not printable
				  	fprintf (stderr,
				    	     "Unknown option character `\\x%x'.\n",
				             optopt);
					cout_usage();
				}
				exit(EXIT_FAILURE);
			default:
				abort ();
		}
    }
}