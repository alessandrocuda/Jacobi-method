#include "barrier.hpp"

barrier::barrier(int n)
    : N_THREADS(n), count(0), generation(0) {}

bool barrier::busywait(std::function<void()> f) {
  int my_gen = generation.load();
  if (count.fetch_add(1) == N_THREADS - 1) {
    f();
    count.store(0);
    generation.fetch_add(1);
    return true;
  } else {
    do {
    } while (my_gen == generation.load());
    return false;
  }
}