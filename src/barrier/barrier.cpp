#include "barrier.hpp"

barrier::barrier(uint64_t count): _count(count){
    _n.store(count);
  }

void barrier::busywait(std::function<void()> f) {
    bool priv_sense = _sense.load();
    if (_n.fetch_sub(1) == 1){
      f();
      _n.store(_count);
      _sense.store(!priv_sense);
    }else{
      while (priv_sense == _sense.load());
    }
}