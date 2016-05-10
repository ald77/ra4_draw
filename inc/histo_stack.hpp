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
#include "TLine.h"
#include "TGraphAsymmErrors.h"

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
    //Probably a better way to do this than using scaled_hist_
    //as mutable storage for formatted plot
    mutable TH1D scaled_hist_;

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
             const std::vector<PlotOpt> &plot_options = {PlotOpt()});
  HistoStack(const HistoStack &) = default;
  HistoStack& operator=(const HistoStack &) = default;
  HistoStack(HistoStack &&) = default;
  HistoStack& operator=(HistoStack &&) = default;
  ~HistoStack() = default;

  void PrintPlot(double luminosity);

  const TH1D & RawHisto(const std::shared_ptr<Process> &process) const;
  TH1D & RawHisto(const std::shared_ptr<Process> &process);

  std::set<std::shared_ptr<Process> > GetProcesses() const;

  std::vector<SingleHist> backgrounds_, signals_, datas_;
  HistoDef definition_;
  std::vector<PlotOpt> plot_options_;

private:
  mutable PlotOpt this_opt_;

  HistoStack() = delete;

  const std::vector<SingleHist> & GetHistoList(const std::shared_ptr<Process> &process) const;
  std::vector<SingleHist> & GetHistoList(const std::shared_ptr<Process> &process);

  const SingleHist & Histo(const std::shared_ptr<Process> &process) const;
  SingleHist & Histo(const std::shared_ptr<Process> &process);

  void RefreshScaledHistos(double luminosity);
  void InitializeHistos() const;
  void MergeOverflow() const;
  void ScaleHistos(double luminosity) const;
  void StackHistos() const;

  void SetRanges() const;

  void ApplyStyles() const;
  void StyleHisto(TH1D &h) const;
  void AdjustFillStyles() const;

  void GetPads(std::unique_ptr<TCanvas> &c,
               std::unique_ptr<TPad> &top,
               std::unique_ptr<TPad> &bottom) const;

  double FixYAxis(std::vector<TH1D> &bottom_plots, TPad *top, TPad *bottom) const;

  std::vector<std::shared_ptr<TLatex> > GetTitleTexts(double luminosity, double left_margin) const;
  TGraphAsymmErrors GetBackgroundError() const;
  std::vector<TLine> GetCutLines() const;
  std::vector<TH1D> GetBottomPlots() const;
  TLine GetBottomHorizontal() const;

  void StripTopPlotLabels() const;

  double GetMaxDraw(double max_bound = std::numeric_limits<double>::infinity()) const;
  double GetMinDraw(double min_bound = 0.) const;

  std::vector<std::shared_ptr<TLegend> > GetLegends(double left_margin);
  void AddEntries(std::vector<std::shared_ptr<TLegend> > &legends,
                  const std::vector<HistoStack::SingleHist> &hists,
                  const std::string &style,
                  std::size_t n_entries,
                  std::size_t &entries_added) const;
  double GetLegendRatio() const;

  double GetYield(std::vector<HistoStack::SingleHist>::const_iterator h) const;
  double GetMean(std::vector<HistoStack::SingleHist>::const_iterator h) const;
};

#endif
