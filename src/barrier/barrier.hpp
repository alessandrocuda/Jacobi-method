#include <atomic>


class barrier {
private:
    const int N_THREADS;
    int counts[2];
    int current;

public:
    barrier(int n);
    std::atomic<int> count;
    std::atomic<int> generation;
    bool busywait();
};

