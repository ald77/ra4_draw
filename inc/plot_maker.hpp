#ifndef H_PLOT_MAKER
#define H_PLOT_MAKER

#include <vector>
#include <set>
#include <memory>
#include <utility>

#include "TCanvas.h"

#include "histo_def.hpp"
#include "process.hpp"
#include "histo_stack.hpp"
#include "plot_opt.hpp"

class PlotMaker{
public:
  PlotMaker();
  PlotMaker(const PlotMaker &) = default;
  PlotMaker& operator=(const PlotMaker &) = default;
  PlotMaker(PlotMaker &&) = default;
  PlotMaker& operator=(PlotMaker &&) = default;
  ~PlotMaker() = default;

  void AddPlot(const std::vector<std::shared_ptr<Process> > & processes,
               const HistoDef &histo_def,
               const PlotOpt &plot_options = PlotOpt());
  
  void MakePlots();

  std::set<std::string> file_extensions_;
private:
  std::vector<HistoStack> stacks_;

  void FillHistograms();
  void FillHistogram(const std::shared_ptr<Process> &process);
  void PrintPlot(HistoStack &stack);
  std::unique_ptr<TLegend> GetLegend(HistoStack &stack);
  std::set<std::shared_ptr<Process> > GetProcesses() const;
  std::vector<std::pair<HistoDef, const TH1D * const> > GetHistos(const std::shared_ptr<Process> &process) const;
  std::vector<std::pair<HistoDef, TH1D * const> > GetHistos(const std::shared_ptr<Process> &process);
};

#endif
