#ifndef H_PROCESS
#define H_PROCESS

#include <memory>
#include <string>

#include "baby.hpp"
#include "named_func.hpp"

class Process : public TAttFill, public TAttLine, public TAttMarker{
public:
  enum class Type{data, background, signal};
  Process(const std::string & name,
          Type type,
          int color,
          std::unique_ptr<Baby> baby,
          const NamedFunc &cut = NamedFunc(1.));
  Process(Process &&) = default;
  Process& operator=(Process &&) = default;
  ~Process() = default;
  
  std::string name_;
  Type type_;
  std::unique_ptr<Baby> baby_;
  NamedFunc cut_;

private:
  Process() = delete;
  Process(const Process &) = delete;
  Process& operator=(const Process &) = delete;
};

#endif
