#ifndef H_HIST1D
#define H_HIST1D

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

#include "core/figure.hpp"
#include "core/process.hpp"
#include "core/axis.hpp"
#include "core/plot_opt.hpp"

class Hist1D final: public Figure{
public:
  class SingleHist1D final: public Figure::FigureComponent{
  public:
    SingleHist1D(const Hist1D &stack,
               const std::shared_ptr<Process> &process,
               const TH1D &hist);
    ~SingleHist1D() = default;

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
    SingleHist1D() = delete;
    SingleHist1D(const SingleHist1D &) = delete;
    SingleHist1D& operator=(const SingleHist1D &) = delete;
    SingleHist1D(SingleHist1D &&) = delete;
    SingleHist1D& operator=(SingleHist1D &&) = delete;

    NamedFunc proc_and_hist_cut_;
    NamedFunc::VectorType cut_vector_, wgt_vector_, val_vector_;
  };

  Hist1D(const Axis &xaxis, const NamedFunc &cut,
             const std::vector<std::shared_ptr<Process> > &processes,
             const std::vector<PlotOpt> &plot_options = {PlotOpt()});
  Hist1D(Hist1D &&) = default;
  Hist1D& operator=(Hist1D &&) = default;
  ~Hist1D() = default;

  void Print(double luminosity,
             const std::string &subdir) final;

  std::set<const Process*> GetProcesses() const final;

  FigureComponent * GetComponent(const Process *process) final;

  std::string Name() const;
  std::string Title() const;

  Hist1D & Weight(const NamedFunc &weight);
  Hist1D & Tag(const std::string &tag);
  Hist1D & LeftLabel(const std::vector<std::string> &label);
  Hist1D & RightLabel(const std::vector<std::string> &label);
  Hist1D & YAxisZoom(const double &yaxis_zoom);
  Hist1D & RatioTitle(const std::string &numerator,
                      const std::string &denominator);

  Axis xaxis_;//!<Specification of content: plotted variable, binning, etc.
  NamedFunc cut_;//!<Event selection
  NamedFunc weight_;//!<Event weight
  std::string tag_;//!<Filename tag to identify plot
  std::vector<std::string> left_label_;//!<Label to plot under the legend, to the left
  std::vector<std::string> right_label_;//!<Label to plot under the legend, to the right
  double yaxis_zoom_;//!<Y-axis zoom
  std::string ratio_numerator_;//!<Label for numerator in ratio plot
  std::string ratio_denominator_;//!<Label for denominator in ratio plot
  std::vector<PlotOpt> plot_options_;//!<Styles with which to draw plot

private:
  std::vector<std::unique_ptr<SingleHist1D> > backgrounds_;//!<Background components of the figure
  std::vector<std::unique_ptr<SingleHist1D> > signals_;//!<Signal components of the figure
  std::vector<std::unique_ptr<SingleHist1D> > datas_;//!<Data components of the figure

  mutable PlotOpt this_opt_;//!<Plot style currently being drawn
  mutable double luminosity_;//!<Luminosity currently being drawn
  mutable double mc_scale_;//!<data/MC normalization
  mutable double mc_scale_error_;//!<data/MC normalization uncertainty
  static TH1D blank_;//<!Blank histogram for creating dummy legend entries

  Hist1D(const Hist1D &) = delete;
  Hist1D& operator=(const Hist1D &) = delete;
  Hist1D() = delete;

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
  std::vector<TLine> GetCutLines(double y_min, double y_max, bool adjust_bottom) const;
  std::vector<TH1D> GetBottomPlots(double &the_min, double &the_max) const;
  TLine GetBottomHorizontal() const;

  void StripTopPlotLabels() const;

  double GetMaxDraw(double max_bound = std::numeric_limits<double>::infinity()) const;
  double GetMinDraw(double min_bound = 0.) const;

  std::vector<std::shared_ptr<TLegend> > GetLegends();
  void AddEntries(std::vector<std::shared_ptr<TLegend> > &legends,
                  const std::vector<std::unique_ptr<SingleHist1D> > &hists,
                  const std::string &style,
                  std::size_t n_entries,
                  std::size_t &entries_added) const;
  double GetLegendRatio() const;

  double GetYield(std::vector<std::unique_ptr<SingleHist1D> >::const_iterator h) const;
  double GetMean(std::vector<std::unique_ptr<SingleHist1D> >::const_iterator h) const;

  void GetTitleSize(double &width, double &height, bool in_pixels) const;

  const std::vector<std::unique_ptr<SingleHist1D> >& GetComponentList(const Process *process);
};

#endif
