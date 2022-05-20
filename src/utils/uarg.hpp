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
					ulong &n, ulong &mode, ulong &nw, unsigned int &seed, int &verbose);

#endif // UARG_H
