#include "process.hpp"

using namespace std;

Process::Process(const string &name,
                 Type type,
                 int color,
                 unique_ptr<Baby> baby,
                 const NamedFunc &cut):
  TAttFill(color, type == Type::background ? 1001 : 0),
  TAttLine(type == Type::signal ? color : 1, 1, 5),
  TAttMarker(color, 20, 1.2),
  name_(name),
  type_(type),
  baby_(move(baby)),
  cut_(cut){
  }
