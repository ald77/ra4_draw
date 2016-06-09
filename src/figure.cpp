#include "figure.hpp"

#include "utilities.hpp"

using namespace std;

Figure::Figure(const vector<shared_ptr<Process> > &/*processes*/):
  backgrounds_(),
  signals_(),
  datas_(){
}

Figure::FigureComponent::FigureComponent(const Figure &figure,
                                         const shared_ptr<Process> &process):
  figure_(figure),
  process_(process){
}

set<shared_ptr<Process> > Figure::GetProcesses() const{
  set<shared_ptr<Process> > processes;
  for(const auto &proc: backgrounds_){
    processes.insert(proc->process_);
  }
  for(const auto &proc: signals_){
    processes.insert(proc->process_);
  }
  for(const auto &proc: datas_){
    processes.insert(proc->process_);
  }
  return processes;
}

Figure::FigureComponent * Figure::GetComponent(const shared_ptr<Process> &process){
  const auto &component_list = GetComponentList(process);
  for(const auto &component: component_list){
    if(component->process_ == process){
      return component.get();
    }
  }
  ERROR("Could not find histogram for process "+process->name_+".");
  return component_list.front().get();
}

vector<unique_ptr<Figure::FigureComponent> >& Figure::GetComponentList(const shared_ptr<Process> &process){
  switch(process->type_){
  case Process::Type::data:
    return datas_;
  case Process::Type::background:
    return backgrounds_;
  case Process::Type::signal:
    return signals_;
  default:
    ERROR("Did not understand process type "+to_string(static_cast<long>(process->type_))+".");
    return backgrounds_;
  }
}
