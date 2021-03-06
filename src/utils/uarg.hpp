#ifndef UARG_H
#define UARG_H

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <ctype.h>
#include "utils.hpp"
#include <string.h>
#include <iostream>
#include <getopt.h>

void 
init_argv(const int argc, char *const argv[],
					uint64_t &n, uint64_t &mode, uint64_t &nw, 
					unsigned int &seed, float &l_range, float &r_range, 
					float &tol, int &verbose);

#endif // UARG_H
