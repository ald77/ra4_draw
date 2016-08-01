#ifndef H_PROCESS
#define H_PROCESS

#include <memory>
#include <string>
#include <map>
#include <set>
#include <mutex>

#include "baby.hpp"
#include "named_func.hpp"
#include "utilities.hpp"

class Process : public TAttFill, public TAttLine, public TAttMarker{
public:
  enum class Type{data, background, signal};

  template<typename BabyType>
  static std::shared_ptr<Process> MakeShared(const std::string &name,
                                             Type type,
                                             int color,
                                             const std::set<std::string> &files,
                                             const NamedFunc &cut = true){
    return std::shared_ptr<Process>(new Process(static_cast<BabyType*>(nullptr),
                                                name, type, color, files, cut));
  }

  std::string name_;
  Type type_;
  NamedFunc cut_;

  std::set<Baby*> Babies() const;

  ~Process();

private:
  template<typename BabyType>
  Process(BabyType * dummy_baby,
          const std::string &name,
          Type type,
          int color,
          const std::set<std::string> &files,
          const NamedFunc &cut);

  Process() = delete;
  Process(const Process &) = delete;
  Process& operator=(const Process &) = delete;
  Process(Process &&) = delete;
  Process& operator=(Process &&) = delete;

  static std::set<std::unique_ptr<Baby> > baby_pool_;
  static std::mutex mutex_;
};

template<typename BabyType>
Process::Process(BabyType * /*dummy_baby*/,
                 const std::string &name,
                 Type type,
                 int color,
                 const std::set<std::string> &files,
                 const NamedFunc &cut):
  TAttFill(color, type == Type::background ? 1001 : 0),
  TAttLine(type == Type::background ? 1 : color,
           1,
           type == Type::background ? 1 :
           type == Type::signal ? 5
           : 3),
  TAttMarker(color, 20, 1.2),
  name_(name),
  type_(type),
  cut_(cut){
  std::lock_guard<std::mutex> lock(mutex_);
  for(const auto &file: files){
    const auto &full_files = Glob(file);
    for(const auto &full_file: full_files){
      bool found = false;
      for(auto &baby_p: baby_pool_){
        auto &baby = *baby_p;
        if(typeid(baby) != typeid(BabyType)) continue;
        const auto &baby_files = baby.FileNames();
        if(baby_files.find(full_file) != baby_files.end()){
          baby.processes_.insert(this);
          found = true;
          break;
        }
      }
      if(!found){
        baby_pool_.emplace(static_cast<Baby*>(new BabyType(std::set<std::string>{full_file},
                                                           std::set<const Process*>{this})));
      }
    }
  }
  }

#endif
