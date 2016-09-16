#ifndef H_PLOT_MAKER
#define H_PLOT_MAKER

#include <vector>
#include <set>
#include <memory>
#include <utility>

#include "core/plot_opt.hpp"
#include "core/figure.hpp"

class Process;

class PlotMaker{
public:
  PlotMaker();
  PlotMaker(const PlotMaker &) = default;
  PlotMaker& operator=(const PlotMaker &) = default;
  PlotMaker(PlotMaker &&) = default;
  PlotMaker& operator=(PlotMaker &&) = default;
  ~PlotMaker() = default;

  template<typename FigureType, typename... Args>
  FigureType & Push(Args&&... args){
    figures_.emplace_back(static_cast<Figure*>(new FigureType(args...)));
    return *static_cast<FigureType*>(figures_.back().get());
  }

  void MakePlots(double luminosity,
                 const std::string &subdir = "");

  const std::vector<std::unique_ptr<Figure> > & Figures() const;
  template<typename FigureType>
  FigureType * GetLast(){
    FigureType *out = nullptr;
    for(auto f = figures_.crbegin(); f != figures_.crend(); ++f){
      if((out = dynamic_cast<FigureType*>(f->get()))) return out;
    }
    return nullptr;
  }
  void Clear();

  bool multithreaded_;
  bool min_print_;

private:
  std::vector<std::unique_ptr<Figure> > figures_;//!<Figures to be produced

  void GetYields();
  long GetYield(Baby *baby_ptr);

  std::set<Baby*> GetBabies() const;
  std::set<const Process *> GetProcesses() const;
  std::set<Figure::FigureComponent*> GetComponents(const Process *process) const;
};

#endif
