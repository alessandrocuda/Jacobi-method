#include "uarg.hpp"

static void
cout_usage()
{
    fputs("Jacobi Test: [OPTIONS]\n", stdout);
    //fputs("\n", stdout);
    fputs("-n   N           define matrix dimension NxN, N >1 \n", stdout);
    fputs("-m   mode        define jacobi implementation mode = [seq, th, ff]\n", stdout);
    fputs("                 seq: sequantion, th: threads, ff: fastflow\n", stdout);
    fputs("-w   nw          define number of workers (nw) in paralell computation\n", stdout);
    fputs("-s   seed        generate a linear system with this seed\n", stdout);
    fputs("-h               print this help and exit\n", stdout);
    fputs("-v               run in verbose mode, for debug purpose\n", stdout);
}

void init_argv(const int argc, char *const argv[], 
					ulong &n, ulong &mode, ulong &nw, 
					unsigned int &seed, float &l_range, float &r_range, 
					float &tol, int &verbose) {
	int 	c; 					// temp arg for getopt, GNU c lib std

	while ((c = getopt(argc, argv, "n:m:w:s:l:r:t:hv:")) != -1){
		switch (c){
			case 'n':
				if (isdigit(*optarg)){
					n = (ulong) atoi(optarg);
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
					nw = (ulong) atoi(optarg);
				}
				break;
			case 's':
				if (isdigit(*optarg)){
					seed = (ulong) atoi(optarg);
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
					verbose = (ulong) atoi(optarg);
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