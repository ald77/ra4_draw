/*! \class PlotMaker

  \brief Organizes efficient production of plots with single loop over each
  process

  \link HistoStack HistoStacks\endlink are added to the PlotMaker using
  PlotMaker::AddPlot(). Once all desired plots have been added, a call to
  PlotMaker::MakePlots() determines the full set of \link Process
  Processes\endlink used by all plots, loops once over each Process to fill all
  histograms using that Process, and then prints the plots.
*/
#include "plot_maker.hpp"

#include <functional>

#include "TLegend.h"

#include "utilities.hpp"
#include "timer.hpp"
#include "thread_pool.hpp"
#include "named_func.hpp"

using namespace std;
using namespace PlotOptTypes;

using ScalarType = NamedFunc::ScalarType;
using VectorType = NamedFunc::VectorType;
using ScalarFunc = NamedFunc::ScalarFunc;
using VectorFunc = NamedFunc::VectorFunc;

/*!\brief Standard constructor
 */
PlotMaker::PlotMaker():
  stacks_(){
}

/*!\brief Add HistoStack to list of plots to be produced at next
  PlotMaker::MakePlots call

  All arguments passed directly to HistoStack constructor. Adds a plot to the
  list of plots to be produced. Does NOT fill or draw the histogram.

  \param[in] histo_def Histogram definition (variable, binning, etc.)

  \param[in] processes Processes to include in plot

  \param[in] plot_options List of styles with which to produce plot
*/
void PlotMaker::AddPlot(const HistoDef &histo_def,
                        const vector<shared_ptr<Process> > &processes,
                        const vector<PlotOpt> &plot_options){
  stacks_.emplace_back(processes, histo_def, plot_options);
}

/*!\brief Prints all added plots with given luminosity

  \param[in] luminosity Integrated luminosity with which to draw plots
*/
void PlotMaker::MakePlots(double luminosity){
  GetYields();

  for(auto &stack: stacks_){
    stack.Print(luminosity);
  }
}

/*!\brief Empties list of plots to be produced at next PlotMaker::MakePlots call
 */
void PlotMaker::Clear(){
  stacks_.clear();
}

void PlotMaker::GetYields(){
  set<shared_ptr<Process> > processes = GetProcesses();
  ThreadPool tp(thread::hardware_concurrency());
  for(const auto &proc: processes){
    tp.Push(bind(&PlotMaker::GetYield, this, ref(proc)));
  }
}

void PlotMaker::GetYield(const std::shared_ptr<Process> &process){
  cout << "Filling histograms for the " << process->name_ << " process..." << endl;

  Baby &baby = *(process->baby_);

  set<Figure::FigureComponent*> figure_components = GetComponents(process);

  long num_entries = baby.GetEntries();
  Timer timer(num_entries, 1.);
  timer.Start();
  for(long entry = 0; entry < num_entries; ++entry){
    timer.Iterate();
    baby.GetEntry(entry);

    for(auto &component: figure_components){
      component->RecordEvent(baby, process->cut_);
    }
  }
}

std::set<std::shared_ptr<Process> > PlotMaker::GetProcesses() const{
  set<shared_ptr<Process> > processes;
  for(const auto &figure: figures_){
    for(const auto &process: figure->GetProcesses()){
      processes.insert(process);
    }
  }
  return processes;
}

std::set<Figure::FigureComponent*> PlotMaker::GetComponents(const std::shared_ptr<Process> &process) const{
  set<Figure::FigureComponent*> figure_components;
  for(auto &figure: figures_){
    auto processes = figure->GetProcesses();
    auto loc = processes.find(process);
    if(loc == processes.end()) continue;
    figure_components.insert(figure->GetComponent(process));
  }
  return figure_components;
}
