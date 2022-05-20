#ifndef UTILS_H
#define UTILS_H

#include <map>
#include <string>
typedef unsigned long ulong;

#define DF_N 10
#define DF_W 1
#define SEQ 1
#define TH 2
#define FF 3

#define SEQUENTIAL "seq"
#define THREADS "th"
#define FASTFLOW "ff"

static const std::map<ulong, const std::string> m = {
                {SEQ, "sequential"}, 
                {TH, "threads"}, 
                {FF, "fastflow"}};


#endif // UTILS_H
