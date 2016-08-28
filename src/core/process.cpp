#include "core/process.hpp"

#include "core/utilities.hpp"

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

Process::~Process(){
  lock_guard<mutex> lock(mutex_);
  auto baby_ptr_iter = baby_pool_.begin();
  while(baby_ptr_iter != baby_pool_.end()){
    auto current_iter = baby_ptr_iter++;
    Baby &baby = *(current_iter->get());
    baby.processes_.erase(this);
    if(baby.processes_.size() == 0){
      baby_pool_.erase(current_iter);
    }
  }
}
