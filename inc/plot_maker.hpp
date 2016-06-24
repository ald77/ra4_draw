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

  template<typename FigureType, typename... Args>
  void Push(Args&&... args){
    figures_.emplace_back(static_cast<Figure*>(new FigureType(args...)));
  }

  void MakePlots(double luminosity);

  const std::vector<std::unique_ptr<Figure> > & Figures() const;
  void Clear();

  bool multithreaded_;

private:
  std::vector<std::unique_ptr<Figure> > figures_;//!<Figures to be produced

  void GetYields();
  long GetYield(const std::shared_ptr<Process> &proc);

  std::set<std::shared_ptr<Process> > GetProcesses() const;
  std::set<Figure::FigureComponent*> GetComponents(const std::shared_ptr<Process> &process) const;
};

#endif
