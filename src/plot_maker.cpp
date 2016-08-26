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
#include <map>
#include <iomanip>  // setw

#include "TLegend.h"

#include "utilities.hpp"
#include "timer.hpp"
#include "thread_pool.hpp"
#include "named_func.hpp"
#include "process.hpp"

using namespace std;
using namespace PlotOptTypes;

using ScalarType = NamedFunc::ScalarType;
using VectorType = NamedFunc::VectorType;
using ScalarFunc = NamedFunc::ScalarFunc;
using VectorFunc = NamedFunc::VectorFunc;

using Clock = chrono::steady_clock;

namespace{
  mutex print_mutex;
}

/*!\brief Standard constructor
 */
PlotMaker::PlotMaker():
  multithreaded_(true),
  min_print_(false),
  figures_(){
}

/*!\brief Prints all added plots with given luminosity

  \param[in] luminosity Integrated luminosity with which to draw plots
*/
void PlotMaker::MakePlots(double luminosity,
                          const string &subdir){
  GetYields();

  for(auto &figure: figures_){
    figure->Print(luminosity, subdir);
  }
}

const vector<unique_ptr<Figure> > & PlotMaker::Figures() const{
  return figures_;
}

/*!\brief Empties list of plots to be produced at next PlotMaker::MakePlots call
 */
void PlotMaker::Clear(){
  figures_.clear();
}

void PlotMaker::GetYields(){
  auto start_time = Clock::now();

  auto babies = GetBabies();
  size_t num_threads = multithreaded_ ? min(babies.size(), static_cast<size_t>(thread::hardware_concurrency())) : 1;
  cout << "Processing " << babies.size() << " babies with " << num_threads << " threads." << endl;

  long num_entries = 0;

  if(multithreaded_ && num_threads>1){
    vector<future<long> > num_entries_future(babies.size());

    ThreadPool tp(num_threads);
    size_t Nbabies = 0;
    for(const auto &baby: babies){
      num_entries_future.at(Nbabies) = tp.Push(bind(&PlotMaker::GetYield, this, ref(baby)));
      ++Nbabies;
    }
    size_t Nfiles=0;
    long printStep=Nbabies/20+1; // Print up to 20 lines of info
    auto start_entries_time = Clock::now();
    for(auto& entries: num_entries_future){
      num_entries += entries.get();
      Nfiles++;
      if(min_print_ && ((Nfiles-1)%printStep==0 || Nfiles==Nbabies)){
	double seconds = chrono::duration<double>(Clock::now()-start_entries_time).count();
	cout<<"Done "<<setw(log10(Nbabies)+1)<<Nfiles<<"/"<<Nbabies<<" files: "<<setw(10)<<AddCommas(num_entries)
	    <<" entries in "<<HoursMinSec(seconds)<<"  ->  "<<setw(5)<<RoundNumber(num_entries/1000.,1,seconds)
	    <<" kHz "<<endl;
      }
    }
  }else{
    for(const auto &baby: babies){
      num_entries += GetYield(ref(baby));
    }
  }
  auto end_time = Clock::now();
  double num_seconds = chrono::duration<double>(end_time-start_time).count();
  if(!min_print_) cout << endl << num_threads << " threads processed "
		       << babies.size() << " babies with "
		       << AddCommas(num_entries) << " events in "
		       << num_seconds << " seconds = "
		       << 0.001*num_entries/num_seconds << " kHz."
		       << endl;
  cout << endl;
}

long PlotMaker::GetYield(Baby *baby_ptr){
  auto start_time = Clock::now();
  Baby &baby = *baby_ptr;
  auto activator = baby.Activate();
  string tag = "";
  if(baby.FileNames().size() == 1){
    tag = Basename(*baby.FileNames().cbegin());
  }else{
    tag = "Baby for processes";
  }
  ostringstream oss;
  oss << " [";
  for(auto proc = baby.processes_.cbegin(); proc != baby.processes_.cend(); ++proc){
    if(proc != baby.processes_.cbegin()) oss << ", ";
    oss << (*proc)->name_;
  }
  oss << "]" << flush;
  tag += oss.str();
  {
    lock_guard<mutex> lock(print_mutex);
    if(!min_print_) cout << "Processing " << tag << endl;
  }

  long num_entries = baby.GetEntries();
  {
    lock_guard<mutex> lock(print_mutex);
    if(!min_print_) cout << tag << " has " << num_entries << " entries." << endl;
  }

  vector<pair<const Process*, set<Figure::FigureComponent*> > > proc_figs(baby.processes_.size());
  size_t iproc = 0;
  for(const auto &proc: baby.processes_){
    proc_figs.at(iproc).first = proc;
    proc_figs.at(iproc).second = GetComponents(proc);
    ++iproc;
  }

  Timer timer(tag, num_entries, 10.);
  for(long entry = 0; entry < num_entries; ++entry){
    if(!min_print_) timer.Iterate();
    baby.GetEntry(entry);

    for(const auto &proc_fig: proc_figs){
      if(proc_fig.first->cut_.IsScalar()){
        if(!proc_fig.first->cut_.GetScalar(baby)) continue;
      }else{
        if(!HavePass(proc_fig.first->cut_.GetVector(baby))) continue;
      }
      for(const auto &component: proc_fig.second){
	lock_guard<mutex> lock(component->mutex_);
        component->RecordEvent(baby);
      }
    }
  }

  auto end_time = Clock::now();
  double num_seconds = chrono::duration<double>(end_time - start_time).count();
  {
    lock_guard<mutex> lock(print_mutex);
    if(!min_print_) cout << "Finished processing " << tag << ". "
			 << num_entries << " events in " << num_seconds << " seconds = "
			 << 0.001*num_entries/num_seconds << " kHz. "<< endl;
  }
  return num_entries;
}

set<Baby*> PlotMaker::GetBabies() const{
  set<Baby*> babies;
  for(auto &proc: GetProcesses()){
    for(const auto &baby: proc->Babies()){
      babies.insert(baby);
    }
  }
  return babies;
}

set<const Process*> PlotMaker::GetProcesses() const{
  set<const Process*> processes;
  for(const auto &figure: figures_){
    for(const auto &process: figure->GetProcesses()){
      processes.insert(process);
    }
  }
  return processes;
}

set<Figure::FigureComponent*> PlotMaker::GetComponents(const Process *process) const{
  set<Figure::FigureComponent*> figure_components;
  for(auto &figure: figures_){
    auto processes = figure->GetProcesses();
    auto loc = processes.find(process);
    if(loc == processes.end()) continue;
    figure_components.insert(figure->GetComponent(process));
  }
  return figure_components;
}
