#include "histo_stack.hpp"

#include <algorithm>

using namespace std;

HistoStack::HistoStack(const std::vector<shared_ptr<Process> > &processes,
                       const HistoDef &definition):
  backgrounds_(),
  signals_(),
  datas_(),
  definition_(definition){
  TH1D empty("",(";"+definition.x_title_+";").c_str(),definition.GetNbins(), &definition.bins_.at(0));
  empty.Sumw2(true);
  for(const auto &process: processes){
    auto ph = make_pair(process, empty);
    ph.second.SetFillColor(process->GetFillColor());
    ph.second.SetFillStyle(process->GetFillStyle());
    ph.second.SetLineColor(process->GetLineColor());
    ph.second.SetLineStyle(process->GetLineStyle());
    ph.second.SetLineWidth(process->GetLineWidth());
    ph.second.SetMarkerColor(process->GetMarkerColor());
    ph.second.SetMarkerStyle(process->GetMarkerStyle());
    ph.second.SetMarkerSize(process->GetMarkerSize());
    switch(process->type_){
    case Process::Type::data:
      datas_.push_back(ph);
      break;
    case Process::Type::background:
      backgrounds_.push_back(ph);
      break;
    case Process::Type::signal:
      signals_.push_back(ph);
      break;
    default:
      break;
    }
  }
}

set<shared_ptr<Process> > HistoStack::GetProcesses() const{
  set<shared_ptr<Process> > processes;
  for(const auto &proc: backgrounds_){
    processes.insert(proc.first);
  }
  for(const auto &proc: signals_){
    processes.insert(proc.first);
  }
  for(const auto &proc: datas_){
    processes.insert(proc.first);
  }
  return processes;
}

const TH1D * HistoStack::GetHisto(const std::shared_ptr<Process> &process) const{
  auto *procs = &backgrounds_;
  switch(process->type_){
  case Process::Type::data:
    procs = &datas_;
    break;
  case Process::Type::background:
    procs = &backgrounds_;
    break;
  case Process::Type::signal:
    procs = &signals_;
    break;
  default:
    procs = nullptr;
    break;
  }
  if(procs == nullptr) return nullptr;
  auto loc = procs->cend();
  for(auto iter = procs->cbegin(); iter != procs->cend(); ++iter){
    if(iter->first == process){
      loc = iter;
      break;
    }
  }
  if(loc == procs->cend()) return nullptr;
  return &(loc->second);
}

TH1D * HistoStack::GetHisto(const std::shared_ptr<Process> &process){
  auto *procs = &backgrounds_;
  switch(process->type_){
  case Process::Type::data:
    procs = &datas_;
    break;
  case Process::Type::background:
    procs = &backgrounds_;
    break;
  case Process::Type::signal:
    procs = &signals_;
    break;
  default:
    procs = nullptr;
    break;
  }
  if(procs == nullptr) return nullptr;
  auto loc = procs->end();
  for(auto iter = procs->begin(); iter != procs->end(); ++iter){
    if(iter->first == process){
      loc = iter;
      break;
    }
  }
  if(loc == procs->end()) return nullptr;
  return &(loc->second);
}
