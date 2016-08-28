#include "core/event_scan.hpp"

#include <iostream>
#include <iomanip>

#include <sys/stat.h>

#include "core/utilities.hpp"

using namespace std;

EventScan::SingleScan::SingleScan(const EventScan &event_scan,
                                  const shared_ptr<Process> &process):
  FigureComponent(event_scan, process),
  out_((event_scan.name_+"_SCAN_"+process->name_+".txt").c_str()),
  full_cut_(event_scan.cut_ && process->cut_),
  cut_vector_(),
  val_vectors_(event_scan.columns_.size()),
  row_(0){
  out_.precision(event_scan.Precision());
}

void EventScan::SingleScan::RecordEvent(const Baby &baby){
  const EventScan &scan = static_cast<const EventScan&>(figure_);
  int w = scan.width_;

  if(full_cut_.IsScalar()){
    if(!full_cut_.GetScalar(baby)) return;
  }else{
    cut_vector_ = full_cut_.GetVector(baby);
  }
  
  size_t max_size = 0;
  for(size_t icol = 0; icol < scan.columns_.size(); ++icol){
    const NamedFunc& col = scan.columns_.at(icol);
    if(col.IsScalar()){
      if(max_size < 1) max_size = 1;
    }else{
      val_vectors_.at(icol) = col.GetVector(baby);
      if(val_vectors_.at(icol).size() > max_size){
	max_size = val_vectors_.at(icol).size();
      }
    }
  }
  if(full_cut_.IsVector() && max_size > cut_vector_.size()){
    max_size = cut_vector_.size();
  }

  if(max_size > 0 && !(row_ & 0x7)){
    out_ << "      Row Instance";
    for(const auto &col: scan.columns_){
      out_ << ' ' << setw(w) << col.Name().substr(0,scan.width_);
    }
    out_.put('\n');
  }

  for(size_t instance = 0; instance < max_size; ++instance){
    out_ << setw(9) << row_ << ' ' << setw(8) << instance;
    for(size_t icol = 0; icol < scan.columns_.size(); ++icol){
      const NamedFunc& col = scan.columns_.at(icol);
      if(col.IsScalar()){
        out_ << ' ' << setw(w) << col.GetScalar(baby);
      }else{
        if(instance < val_vectors_.at(icol).size()){
          out_ << ' ' << setw(w) << val_vectors_.at(icol).at(instance);
        }else{
	  out_ << ' ' << setw(w) << ' ';
	}
      }
    }
    out_.put('\n');
  }

  if(max_size > 0) ++row_;
}

void EventScan::SingleScan::Precision(unsigned precision){
  out_.precision(precision);
}

EventScan::EventScan(const string &name,
                     const NamedFunc &cut,
                     const vector<NamedFunc> &columns,
                     const vector<shared_ptr<Process> > &processes,
		     unsigned precision):
  name_(name),
  cut_(cut),
  columns_(columns),
  scans_(),
  precision_(precision),
  width_(precision+6){
  for(const auto& proc: processes){
    scans_.emplace_back(new SingleScan(*this, proc));
  }
}

void EventScan::Print(double /*luminosity*/,
                      const std::string & /*subdir*/){
  for(const auto &scan: scans_){
    cout
      << "Wrote scan to "
      << name_
      << "_SCAN_"
      << scan->process_->name_ << ".txt"
      << endl;
  }
}

set<const Process*> EventScan::GetProcesses() const{
  set<const Process *> processes;
  for(const auto &scan: scans_){
    processes.insert(scan->process_.get());
  }
  return processes;
}

Figure::FigureComponent * EventScan::GetComponent(const Process *process){
  for(const auto &scan: scans_){
    if(scan->process_.get() == process) return scan.get();
  }
  return nullptr;
}

unsigned EventScan::Precision() const{
  return precision_;
}

EventScan & EventScan::Precision(unsigned precision){
  precision_ = precision;
  width_ = precision + 6;
  for(auto &scan: scans_){
    scan->Precision(precision);
  }
  return *this;
}
