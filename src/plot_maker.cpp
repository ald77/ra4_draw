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
#include <mutex>
#include <chrono>

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

using Clock = chrono::high_resolution_clock;

namespace{
  mutex print_mutex;
}

/*!\brief Standard constructor
 */
PlotMaker::PlotMaker():
  multithreaded_(true),
  figures_(){
}

/*!\brief Prints all added plots with given luminosity

  \param[in] luminosity Integrated luminosity with which to draw plots
*/
void PlotMaker::MakePlots(double luminosity){
  GetYields();

  for(auto &figure: figures_){
    figure->Print(luminosity);
  }
}

const std::vector<std::unique_ptr<Figure> > & PlotMaker::Figures() const{
  return figures_;
}

/*!\brief Empties list of plots to be produced at next PlotMaker::MakePlots call
 */
void PlotMaker::Clear(){
  figures_.clear();
}

void PlotMaker::GetYields(){
  auto start_time = Clock::now();

  set<shared_ptr<Process> > processes = GetProcesses();
  size_t num_threads = multithreaded_ ? min(processes.size(), static_cast<size_t>(thread::hardware_concurrency())) : 1;
  cout << "Processing " << processes.size() << " samples with " << num_threads << " threads." << endl;

  long num_entries = 0;

  if(multithreaded_ && num_threads>1){
    vector<future<long> > num_entries_future(processes.size());

    ThreadPool tp(num_threads);
    size_t i = 0;
    for(const auto &proc: processes){
      num_entries_future.at(i) = tp.Push(bind(&PlotMaker::GetYield, this, ref(proc)));
      ++i;
    }
    for(auto& entries: num_entries_future){
      num_entries += entries.get();
    }
  }else{
    for(const auto &proc: processes){
      num_entries += GetYield(ref(proc));
    }
  }

  auto end_time = Clock::now();
  double num_seconds = chrono::duration<double>(end_time-start_time).count();
  cout << num_threads << " threads finished "
       << processes.size() << " samples with "
       << num_entries << " events in "
       << num_seconds << " seconds = "
       << 0.001*num_entries/num_seconds << " kHz."
       << endl;
}

long PlotMaker::GetYield(const std::shared_ptr<Process> &process){
  auto start_time = Clock::now();
  {
    lock_guard<mutex> lock(print_mutex);
    cout << "Starting to process " << process->name_ << "." << endl;
  }

  Baby &baby = *(process->baby_);

  long num_entries = baby.GetEntries();
  {
    lock_guard<mutex> lock(print_mutex);
    cout << process->name_ << " has " << num_entries << " entries and cut "<<process->cut_.PlainName() << endl;
  }

  set<Figure::FigureComponent*> figure_components = GetComponents(process);

  Timer timer(num_entries, 10.);
  timer.Start();
  for(long entry = 0; entry < num_entries; ++entry){
    timer.Iterate();
    baby.GetEntry(entry);

    if(process->cut_.IsScalar()){
      if(!process->cut_.GetScalar(baby)) continue;
    }else{
      if(!HavePass(process->cut_.GetVector(baby))) continue;
    }

    for(auto &component: figure_components){
      component->RecordEvent(baby);
    }
  }

  auto end_time = Clock::now();
  double num_seconds = chrono::duration<double>(end_time - start_time).count();
  {
    lock_guard<mutex> lock(print_mutex);
    cout << "Finished processing " << process->name_ << ". "
	 << num_entries << " events in " << num_seconds << " seconds = "
	 << 0.001*num_entries/num_seconds << " kHz. "<< endl;
  }
  return num_entries;
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
