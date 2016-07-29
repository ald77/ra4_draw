#include "process.hpp"

#include "utilities.hpp"

using namespace std;

set<unique_ptr<Baby> > Process::baby_pool_{};
mutex Process::mutex_{};

set<Baby*> Process::Babies() const{
  lock_guard<mutex> lock(mutex_);
  set<Baby*> babies;
  for(const auto &baby: baby_pool_){
    const auto &procs = baby->processes_;
    if(procs.find(this) != procs.end()){
      babies.insert(baby.get());
    }
  }
  return babies;
}
