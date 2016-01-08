#ifndef H_HISTO_STACK
#define H_HISTO_STACK

#include <vector>
#include <utility>
#include <memory>
#include <set>

#include "TH1D.h"

#include "process.hpp"
#include "histo_def.hpp"

class HistoStack{
public:
  HistoStack(const std::vector<std::shared_ptr<Process> > &processes,
             const HistoDef &definition);
  HistoStack(const HistoStack &) = default;
  HistoStack& operator=(const HistoStack &) = default;
  HistoStack(HistoStack &&) = default;
  HistoStack& operator=(HistoStack &&) = default;
  ~HistoStack() = default;

  std::set<std::shared_ptr<Process> > GetProcesses() const;
  const TH1D * GetHisto(const std::shared_ptr<Process> &process) const;
  TH1D * GetHisto(const std::shared_ptr<Process> &process);

  std::vector<std::pair<std::shared_ptr<Process>, TH1D> > backgrounds_, signals_, datas_;
  HistoDef definition_;

private:
  HistoStack() = delete;
};

#endif
