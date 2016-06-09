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

  void AddPlot(const HistoDef &histo_def,
               const std::vector<std::shared_ptr<Process> > & processes,
               const std::vector<PlotOpt> &plot_options = {PlotOpt()});

  void MakePlots(double luminosity);

  void Clear();

private:
  std::vector<HistoStack> stacks_;//!<Plots to be produced

  std::vector<std::shared_ptr<Figure> > figures_;//!<Figures to be produced

  void GetYields();
  void GetYield(const std::shared_ptr<Process> &proc);

  std::set<std::shared_ptr<Process> > GetProcesses() const;
  std::set<Figure::FigureComponent*> GetComponents(const std::shared_ptr<Process> &process) const;
};

#endif
