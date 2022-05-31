#include <atomic>
#include <functional>


class barrier {
private:
    const int N_THREADS;
    int counts[2];

public:
    barrier(int n);
    std::atomic<int> count;
    std::atomic<int> generation;
    bool busywait(std::function<void()> f);
};