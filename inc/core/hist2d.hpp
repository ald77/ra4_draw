#ifndef H_HIST2D
#define H_HIST2D

#include "TH2D.h"
#include "TGraph.h"
#include "TLine.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "core/figure.hpp"
#include "core/axis.hpp"
#include "core/plot_opt.hpp"
#include "core/clusterizer.hpp"

class Hist2D: public Figure{
public:
  class SingleHist2D: public Figure::FigureComponent{
  public:
    SingleHist2D(const Hist2D &figure,
                 const std::shared_ptr<Process> &process,
                 const TH2D &hist_template);
    ~SingleHist2D() = default;

    Clustering::Clusterizer clusterizer_;

    void RecordEvent(const Baby &baby);

  private:
    SingleHist2D() = delete;
    SingleHist2D(const SingleHist2D &) = delete;
    SingleHist2D& operator=(const SingleHist2D &) = delete;
    SingleHist2D(SingleHist2D &&) = delete;
    SingleHist2D& operator=(SingleHist2D &&) = delete;

    NamedFunc proc_and_hist_cut_;
    NamedFunc::VectorType cut_vector_, wgt_vector_, xval_vector_, yval_vector_;
  };

  Hist2D(const Axis &xaxis, const Axis &yaxis, const NamedFunc &cut,
         const std::vector<std::shared_ptr<Process> > &processes,
         const std::vector<PlotOpt> &plot_options = {PlotOpt()});
  Hist2D(Hist2D &&) = default;
  Hist2D& operator=(Hist2D &&) = default;
  ~Hist2D() = default;

  void Print(double luminosity,
             const std::string &subdir) override;

  std::set<const Process*> GetProcesses() const override;

  FigureComponent * GetComponent(const Process *process) override;

  std::string Name() const;

  Hist2D & Weight(const NamedFunc &weight);
  Hist2D & Tag(const std::string &tag);

  Axis xaxis_, yaxis_;
  NamedFunc cut_, weight_;
  std::string tag_;
  std::vector<PlotOpt> plot_options_;

private:
  std::vector<std::unique_ptr<SingleHist2D> > backgrounds_;
  std::vector<std::unique_ptr<SingleHist2D> > signals_;
  std::vector<std::unique_ptr<SingleHist2D> > datas_;

  mutable PlotOpt this_opt_;
  mutable double luminosity_;
  static TH2D blank_;

  Hist2D(const Hist2D &) = delete;
  Hist2D& operator=(const Hist2D &) = delete;
  Hist2D() = delete;

  void MakeOnePlot(const std::string &subdir);
  TH2D GetBkgHist(bool bkg_is_hist) const;
  std::vector<TGraph> GetGraphs(const std::vector<std::unique_ptr<SingleHist2D> > &components,
				bool lumi_weighted) const;
  std::vector<TLine> GetLines() const;
  std::vector<std::shared_ptr<TLatex> > GetLabels(bool bkg_is_hist) const;
  void AddEntry(TLegend &l, const SingleHist2D &h, const TGraph &g) const;

  const std::vector<std::unique_ptr<SingleHist2D> >& GetComponentList(const Process *process);
};

void SetStyle(TPaveText &pt);

#endif
