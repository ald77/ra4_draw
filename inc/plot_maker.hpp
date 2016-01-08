#ifndef H_PLOT_MAKER
#define H_PLOT_MAKER

#include <vector>
#include <set>
#include <memory>
#include <utility>

#include "histo_def.hpp"
#include "Process.hpp"
#include "histo_stack.hpp"

class PlotMaker{
public:
  enum class UpperType{lumiNorm, dataNorm, unstacked, shapeCompare};
  enum class LowerType{none, ratio, diff, residual};
  
  PlotMaker();
  PlotMaker(const PlotMaker &) = default;
  PlotMaker& operator=(const PlotMaker &) = default;
  PlotMaker(PlotMaker &&) = default;
  PlotMaker& operator=(PlotMaker &&) = default;
  ~PlotMaker() = default;

  void AddPlot(const std::vector<std::shared_ptr<Process> > & processes,
               const HistoDef &histo_def);
  
  void MakePlots();

  std::set<std::string> file_extensions_;
  UpperType upper_type_;
  LowerType lower_type_;
  bool do_log_upper_, do_log_lower_;
private:
  std::vector<HistoStack> stacks_;

  void FillHistograms();
  void PrintPlots();
  std::set<std::shared_ptr<Process> > GetProcesses() const;
  std::vector<std::pair<HistoDef, const TH1D * const> > GetHistos(const std::shared_ptr<Process> &process) const;
  std::vector<std::pair<HistoDef, TH1D * const> > GetHistos(const std::shared_ptr<Process> &process);
};

#endif
