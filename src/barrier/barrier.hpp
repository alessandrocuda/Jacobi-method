#include <atomic>
#include <stdint.h>
#include <thread>
#include <stdexcept>
#include <functional>

class barrier {
  uint64_t _count;
  std::atomic_bool _sense;
  std::atomic<uint64_t> _n;

public:
    barrier(uint64_t count);
    void busywait(std::function<void()> f);
};