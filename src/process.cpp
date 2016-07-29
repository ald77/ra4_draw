#include "process.hpp"

#include "utilities.hpp"

using namespace std;

std::set<std::unique_ptr<Baby> > Process::baby_pool_{};

set<Baby*> Process::Babies() const{
  set<Baby*> babies;
  for(const auto &baby: baby_pool_){
    const auto &procs = baby->processes_;
    if(procs.find(this) != procs.end()){
      babies.insert(baby.get());
    }
  }
  return babies;
}

const std::set<std::unique_ptr<Baby> > & Process::BabyPool(){
  return baby_pool_;
}
