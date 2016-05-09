#ifndef H_PLOT_MAKER
#define H_PLOT_MAKER

#include <vector>
#include <set>
#include <memory>
#include <utility>

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
               const std::vector<PlotOpt> &plot_options = {PlotOpt()});

  void MakePlots(double luminosity);

  void Clear();

private:
  std::vector<HistoStack> stacks_;

  void FillHistograms();
  void FillHistogram(const std::shared_ptr<Process> &process);
  std::set<std::shared_ptr<Process> > GetProcesses() const;
  std::vector<std::pair<HistoDef, const TH1D * const> > GetHistos(const std::shared_ptr<Process> &process) const;
  std::vector<std::pair<HistoDef, TH1D * const> > GetHistos(const std::shared_ptr<Process> &process);
};

#endif
