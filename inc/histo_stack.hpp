#ifndef H_HISTO_STACK
#define H_HISTO_STACK

#include <vector>
#include <utility>
#include <memory>
#include <set>
#include <limits>

#include "TH1D.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"

#include "process.hpp"
#include "histo_def.hpp"
#include "plot_opt.hpp"

class HistoStack{
public:
  class SingleHist{
  public:
    SingleHist(const std::shared_ptr<Process> &process,
               const TH1D &hist);
    SingleHist(const SingleHist &) = default;
    SingleHist& operator=(const SingleHist &) = default;
    SingleHist(SingleHist &&) = default;
    SingleHist& operator=(SingleHist &&) = default;
    ~SingleHist() = default;

    std::shared_ptr<Process> process_;
    TH1D raw_hist_;
    TH1D scaled_hist_;

    double GetMax(double max_bound = std::numeric_limits<double>::infinity(),
                  bool include_error_bar = false,
                  bool include_overflow = false) const;
    double GetMin(double max_bound = 0.,
                  bool include_error_bar = false,
                  bool include_overflow = false) const;

  private:
    SingleHist() = delete;
  };

  HistoStack(const std::vector<std::shared_ptr<Process> > &processes,
             const HistoDef &definition,
             const PlotOpt &plot_options = PlotOpt());
  HistoStack(const HistoStack &) = default;
  HistoStack& operator=(const HistoStack &) = default;
  HistoStack(HistoStack &&) = default;
  HistoStack& operator=(HistoStack &&) = default;
  ~HistoStack() = default;

  void StripTopPlotLabels();
  void GetPads(std::unique_ptr<TCanvas> &c,
               std::unique_ptr<TPad> &top,
               std::unique_ptr<TPad> &bottom) const;
  std::vector<std::shared_ptr<TLatex> > GetTitleTexts(double luminosity) const;
  void PrintPlot(double luminosity);

  const TH1D & RawHisto(const std::shared_ptr<Process> &process) const;
  TH1D & RawHisto(const std::shared_ptr<Process> &process);

  const TH1D & ScaledHisto(const std::shared_ptr<Process> &process) const;

  HistoStack & SetPlotOptions(const PlotOpt &plot_opt);
  const PlotOpt & GetPlotOptions() const;

  void RefreshScaledHistos(double luminosity);

  std::vector<TH1D> GetBottomPlots() const;

  std::set<std::shared_ptr<Process> > GetProcesses() const;
  std::unique_ptr<TLegend> GetLegend();

  std::vector<SingleHist> backgrounds_, signals_, datas_;
  HistoDef definition_;
  PlotOpt plot_options_;

private:
  HistoStack() = delete;

  void StackHistos(double luminosity);
  void MergeOverflow();
  void SetRanges();

  const std::vector<SingleHist> & GetHistoList(const std::shared_ptr<Process> &process) const;
  std::vector<SingleHist> & GetHistoList(const std::shared_ptr<Process> &process);

  const SingleHist & Histo(const std::shared_ptr<Process> &process) const;
  SingleHist & Histo(const std::shared_ptr<Process> &process);

  double GetMaxDraw(double max_bound = std::numeric_limits<double>::infinity()) const;
  double GetMinDraw(double min_bound = 0.) const;
};

#endif
