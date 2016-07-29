#ifndef H_EVENT_SCAN
#define H_EVENT_SCAN

#include <memory>
#include <vector>
#include <string>
#include <fstream>

#include "figure.hpp"
#include "process.hpp"

class EventScan final : public Figure{
 public:
  class SingleScan final : public Figure::FigureComponent{
 public:
   SingleScan(const EventScan &event_scan,
              const std::shared_ptr<Process> &process);
   ~SingleScan() = default;

   void RecordEvent(const Baby &baby) final;

   void Precision(unsigned precision);

 private:
   SingleScan() = delete;
   SingleScan(const SingleScan &) = delete;
   SingleScan& operator=(const SingleScan &) = delete;
   SingleScan(SingleScan &&) = delete;
   SingleScan& operator=(SingleScan &&) = delete;

   std::ofstream out_;//!<File to which results are printed
   NamedFunc full_cut_;//!<Cached scan&&process cut
   NamedFunc::VectorType cut_vector_;//!<Cut results (to avoid creating new vector each event)
   std::vector<NamedFunc::VectorType> val_vectors_;//!<Values for each column (to avoid creating new vectors each event)
   std::size_t row_;
 };

 EventScan(const std::string &name,
           const NamedFunc &cut,
           const std::vector<NamedFunc> &columns,
           const std::vector<std::shared_ptr<Process> > &processes,
	   unsigned precision = 10);
 EventScan(EventScan &&) = default;
 EventScan& operator=(EventScan &&) = default;
 ~EventScan() = default;

 void Print(double luminosity,
            const std::string &subdir) final;

 std::set<const Process*> GetProcesses() const final;

 FigureComponent * GetComponent(const Process *process) final;

 unsigned Precision() const;
 EventScan & Precision(unsigned precision);

 std::string name_;//!<Name of scan for saving to file
 NamedFunc cut_;//!<Cut restricting printed events/objects
 std::vector<NamedFunc> columns_;//!<Variables to print

 private:
 std::vector<std::unique_ptr<SingleScan> > scans_;//!<One scan for each process
 unsigned precision_;//!<Decimal places to print
 unsigned width_;//!<Width of column in characters. Determined from precision
 
 EventScan(const EventScan &) = delete;
 EventScan& operator=(const EventScan &) = delete;
 EventScan() = delete;
};

#endif
