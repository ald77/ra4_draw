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

#include "figure.hpp"
#include "process.hpp"
#include "histo_def.hpp"
#include "plot_opt.hpp"

class HistoStack final: public Figure{
public:
  class SingleHist final: public Figure::FigureComponent{
  public:
    SingleHist(const HistoStack &stack,
               const std::shared_ptr<Process> &process,
               const TH1D &hist);
    SingleHist(const SingleHist &) = default;
    SingleHist& operator=(const SingleHist &) = default;
    SingleHist(SingleHist &&) = default;
    SingleHist& operator=(SingleHist &&) = default;
    ~SingleHist() = default;

    TH1D raw_hist_;//!<Histogram storing distribution before stacking and luminosity weighting
    mutable TH1D scaled_hist_;//!<Kludge. Mutable storage of scaled and stacked histogram

    void RecordEvent(const Baby &baby) final;

    double GetMax(double max_bound = std::numeric_limits<double>::infinity(),
                  bool include_error_bar = false,
                  bool include_overflow = false) const;
    double GetMin(double max_bound = 0.,
                  bool include_error_bar = false,
                  bool include_overflow = false) const;

  private:
    SingleHist() = delete;

    NamedFunc proc_and_hist_cut_;
    NamedFunc::VectorType cut_vector_, wgt_vector_, val_vector_;
  };

  HistoStack(const HistoDef &definition,
             const std::vector<std::shared_ptr<Process> > &processes,
             const std::vector<PlotOpt> &plot_options = {PlotOpt()});
  HistoStack(HistoStack &&) = default;
  HistoStack& operator=(HistoStack &&) = default;
  ~HistoStack() = default;

  void Print(double luminosity) final;

  std::set<std::shared_ptr<Process> > GetProcesses() const final;

  FigureComponent * GetComponent(const std::shared_ptr<Process> &process) final;

  HistoDef definition_;//!<Specification of content: plotted variable, binning, etc.
  std::vector<PlotOpt> plot_options_;//!<Styles with which to draw plot

private:
  std::vector<std::unique_ptr<SingleHist> > backgrounds_;//!<Background components of the figure
  std::vector<std::unique_ptr<SingleHist> > signals_;//!<Signal components of the figure
  std::vector<std::unique_ptr<SingleHist> > datas_;//!<Data components of the figure

  mutable PlotOpt this_opt_;//!<Plot style currently being drawn
  mutable double luminosity_;//!<Luminosity currently being drawn
  mutable double mc_scale_;//!<data/MC normalization
  mutable double mc_scale_error_;//!<data/MC normalization uncertainty
  static TH1D blank_;//<!Blank histogram for creating dummy legend entries

  HistoStack(const HistoStack &) = delete;
  HistoStack& operator=(const HistoStack &) = delete;
  HistoStack() = delete;

  const std::vector<SingleHist> & GetHistoList(const std::shared_ptr<Process> &process) const;
  std::vector<SingleHist> & GetHistoList(const std::shared_ptr<Process> &process);

  const SingleHist & Histo(const std::shared_ptr<Process> &process) const;
  SingleHist & Histo(const std::shared_ptr<Process> &process);

  void RefreshScaledHistos();
  void InitializeHistos() const;
  void MergeOverflow() const;
  void ScaleHistos() const;
  void StackHistos() const;
  void NormalizeHistos() const;

  void SetRanges() const;

  void ApplyStyles() const;
  void StyleHisto(TH1D &h) const;
  void AdjustFillStyles() const;

  void GetPads(std::unique_ptr<TCanvas> &c,
               std::unique_ptr<TPad> &top,
               std::unique_ptr<TPad> &bottom) const;

  void FixYAxis(std::vector<TH1D> &bottom_plots) const;

  std::vector<std::shared_ptr<TLatex> > GetTitleTexts() const;
  TGraphAsymmErrors GetBackgroundError() const;
  std::vector<TLine> GetCutLines(double y_min, double y_max) const;
  std::vector<TH1D> GetBottomPlots(double &the_min, double &the_max) const;
  TLine GetBottomHorizontal() const;

  void StripTopPlotLabels() const;

  double GetMaxDraw(double max_bound = std::numeric_limits<double>::infinity()) const;
  double GetMinDraw(double min_bound = 0.) const;

  std::vector<std::shared_ptr<TLegend> > GetLegends();
  void AddEntries(std::vector<std::shared_ptr<TLegend> > &legends,
                  const std::vector<std::unique_ptr<SingleHist> > &hists,
                  const std::string &style,
                  std::size_t n_entries,
                  std::size_t &entries_added) const;
  double GetLegendRatio() const;

  double GetYield(std::vector<std::unique_ptr<SingleHist> >::const_iterator h) const;
  double GetMean(std::vector<std::unique_ptr<SingleHist> >::const_iterator h) const;

  void GetTitleSize(double &width, double &height, bool in_pixels) const;

  const std::vector<std::unique_ptr<SingleHist> >& GetComponentList(const std::shared_ptr<Process> &process);
};

#endif
