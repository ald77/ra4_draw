/*! \class HistoStack

  \brief A full 1D plot with stacked/overlayed histograms

  HistoStack contains all the information necessary to produce a single 1D plot
  containing a combination of background MC, signal MC, and data histograms. The
  content and style of the plot are maintained (mostly) independently, so that
  once the histograms have been filled with data, redrawing with multiple styles
  has minimal overhead.

  To produce a plot, the component histograms in HistoStack::backgrounds_,
  HistoStack::signals_, and HistoStack::datas_ must be filled (e.g., by
  PlotMaker). Once the data is ready, a call to HistoStack::Print() will
  generate the formatted plot for each style contained in
  HistoStack::plot_options_.

*/

/*! \class HistoStack::SingleHist

  \brief Container for a TH1D associated with a single Process

  HistoStack::SingleHist is mostly a "dumb" container used by HistoStack for
  convenience. It contains a pointer to a single Process and a TH1D.

  Until I have a more elegant solution, it also contains a second TH1D which is
  (ab)used by HistoStack to draw the stacked and luminosity scaled histogram
  without disturbing the data in the main TH1D.
*/

#include "core/histo_stack.hpp"

#include <cmath>

#include <algorithm>
#include <sstream>

#include <sys/stat.h>

#include "TROOT.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"
#include "TLegendEntry.h"

#include "core/utilities.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace{
  class Counter{
  public:
    Counter():
      count_(0){
    }

    string operator()(){
      return to_string(count_++);
    }
  private:
    unsigned long count_;
  } counter;

  /*!\brief Draws all histograms to current canvas, updating draw_opt to contain
    "same" as needed

    \param[in] hists List of histograms to draw

    \param[in,out] draw_opt Option string used by TH1D to draw
    histogram. Changed to "hist same" after first histogram is drawn
  */
  void DrawAll(const vector<unique_ptr<HistoStack::SingleHist> > &hists,
               string &draw_opt, bool reversed = false){
    if(!reversed){
      for(auto &hist: hists){
        hist->scaled_hist_.Draw(draw_opt.c_str());
        if(!Contains(draw_opt, "same")){
          draw_opt = draw_opt + " same";
        }
      }
    }else{
      for(auto h = hists.crbegin(); h != hists.crend(); ++h){
        auto &hist = *h;
        hist->scaled_hist_.Draw(draw_opt.c_str());
        if(!Contains(draw_opt, "same")){
          draw_opt = draw_opt + " same";
        }
      }
    }
  }

  /*!\brief Erases x-axis title and labels from plot

    \param[in,out] h Histogram to modify
  */
  void StripXLabels(TH1D &h){
    h.GetXaxis()->SetTitle("");
    h.SetLabelSize(0., "x");
    h.SetTitleSize(0., "x");
  }

  /*!\brief Determines which legend column to put an entry in

    \param[in] entry Number of entries already in legend

    \param[in] n_entries Total entries that will go in legend

    \param[in] n_columns Number of columns in legend

    \return Which column to put entry in. 0 is leftmost, n_columns-1 is
    rightmost.
  */
  size_t GetLegendIndex(size_t entry, size_t n_entries, size_t n_columns){
    size_t entries_per_column = n_entries / n_columns;
    size_t cols_with_extra_entry = n_entries % n_columns;
    size_t this_col = -1;
    size_t this_col_end = 0;
    while(this_col_end <= entry){
      ++this_col;
      this_col_end += entries_per_column;
      if(this_col < cols_with_extra_entry){
        ++this_col_end;
      }
    }
    return this_col;
  }
}

TH1D HistoStack::blank_ = TH1D();

/*!\brief Standard constructor

  \param[in] process Process used to fill histogram

  \param[in] hist A fully styled, and typically unfilled histogram to start from
*/
HistoStack::SingleHist::SingleHist(const HistoStack &figure,
                                   const std::shared_ptr<Process> &process,
                                   const TH1D &hist):
  FigureComponent(figure, process),
  raw_hist_(hist),
  scaled_hist_(),
  proc_and_hist_cut_(figure.definition_.cut_ && process->cut_),
  cut_vector_(),
  wgt_vector_(),
  val_vector_(){
  raw_hist_.Sumw2();
  scaled_hist_.Sumw2();
}

void HistoStack::SingleHist::RecordEvent(const Baby &baby){
  const HistoStack& stack = static_cast<const HistoStack&>(figure_);
  size_t min_vec_size;
  bool have_vec = false;

  const NamedFunc &cut = proc_and_hist_cut_;
  if(cut.IsScalar()){
    if(!cut.GetScalar(baby)) return;
  }else{
    cut_vector_ = cut.GetVector(baby);
    if(!HavePass(cut_vector_)) return;
    have_vec = true;
    min_vec_size = cut_vector_.size();
  }

  const NamedFunc &wgt = stack.definition_.weight_;
  NamedFunc::ScalarType wgt_scalar = 0.;
  if(wgt.IsScalar()){
    wgt_scalar = wgt.GetScalar(baby);
  }else{
    wgt_vector_ = wgt.GetVector(baby);
    if(!have_vec || wgt_vector_.size() < min_vec_size){
      have_vec = true;
      min_vec_size = wgt_vector_.size();
    }
  }

  const NamedFunc &val = stack.definition_.var_;
  NamedFunc::ScalarType val_scalar = 0.;
  if(val.IsScalar()){
    val_scalar = val.GetScalar(baby);
  }else{
    val_vector_ = val.GetVector(baby);
    if(!have_vec || val_vector_.size() < min_vec_size){
      have_vec = true;
      min_vec_size = val_vector_.size();
    }
  }

  if(!have_vec){
    raw_hist_.Fill(val_scalar, wgt_scalar);
  }else{
    for(size_t i = 0; i < min_vec_size; ++i){
      if(cut.IsVector() && !cut_vector_.at(i)) continue;
      raw_hist_.Fill(val.IsScalar() ? val_scalar : val_vector_.at(i),
                     wgt.IsScalar() ? wgt_scalar : wgt_vector_.at(i));
    }
  }
}

/*! Get the maximum of the histogram

  \param[in] max_bound Returns the highest bin content c satisfying
  c<max_bound. Usually infinity.

  \param[in] include_error_bar If true, use bin content+error instead of just
  content

  \param[in] include_overflow If true, also check height of underflow and
  overflow bins

  \return Histogram maximum
*/
double HistoStack::SingleHist::GetMax(double max_bound,
                                      bool include_error_bar,
                                      bool include_overflow) const{
  int start_bin = include_overflow ? 0 : 1;
  int end_bin = include_overflow ? (scaled_hist_.GetNbinsX()+1) : scaled_hist_.GetNbinsX();
  double the_max = -numeric_limits<double>::infinity();
  for(int bin = start_bin; bin <= end_bin; ++bin){
    double content = scaled_hist_.GetBinContent(bin);
    if(include_error_bar){
      content += scaled_hist_.GetBinErrorUp(bin);
    }
    if(content > the_max && content < max_bound){
      the_max = content;
    }
  }
  return the_max;
}

/*! Get the minimum of the histogram

  \param[in] min_bound Returns the lowest bin content c satisfying
  c>min_bound. Usually zero.

  \param[in] include_error_bar If true, use bin content-error instead of just
  content

  \param[in] include_overflow If true, also check height of underflow and
  overflow bins

  \return Histogram minimum
*/
double HistoStack::SingleHist::GetMin(double min_bound,
                                      bool include_error_bar,
                                      bool include_overflow) const{
  int start_bin = include_overflow ? 0 : 1;
  int end_bin = include_overflow ? scaled_hist_.GetNbinsX() : (scaled_hist_.GetNbinsX()+1);
  double the_min = numeric_limits<double>::infinity();
  for(int bin = start_bin; bin <= end_bin; ++bin){
    double content = scaled_hist_.GetBinContent(bin);
    if(include_error_bar){
      content -= fabs(scaled_hist_.GetBinErrorLow(bin));
    }
    if(content < the_min && content > min_bound){
      the_min = content;
    }
  }
  return the_min;
}

/*! \brief Standard constructor

  \param[in] processes List of process for the component histograms

  \param[in] definition Specification of contents (plotted variable, binning,
  etc.)

  \param[in] plot_options Styles with which to draw plot
*/
HistoStack::HistoStack(const HistoDef &definition,
                       const std::vector<std::shared_ptr<Process> > &processes,
                       const std::vector<PlotOpt> &plot_options):
  Figure(),
  definition_(definition),
  plot_options_(plot_options),
  this_opt_(PlotOpt()),
  luminosity_(),
  mc_scale_(),
  mc_scale_error_(){
  if(plot_options_.size() > 0) this_opt_ = plot_options_.front();

  string x_title = definition.x_title_;
  if(definition.units_ != "") x_title += " ["+definition.units_+"]";

  TH1D empty("", (";"+x_title+";").c_str(), definition.Nbins(), &definition.Bins().at(0));
  empty.SetStats(false);
  empty.Sumw2(true);
  for(const auto &process: processes){
    unique_ptr<SingleHist> hist(new SingleHist(*this, process, empty));
    hist->raw_hist_.SetFillColor(process->GetFillColor());
    hist->raw_hist_.SetFillStyle(process->GetFillStyle());
    hist->raw_hist_.SetLineColor(process->GetLineColor());
    hist->raw_hist_.SetLineStyle(process->GetLineStyle());
    hist->raw_hist_.SetLineWidth(process->GetLineWidth());
    hist->raw_hist_.SetMarkerColor(process->GetMarkerColor());
    hist->raw_hist_.SetMarkerStyle(process->GetMarkerStyle());
    hist->raw_hist_.SetMarkerSize(process->GetMarkerSize());

    switch(process->type_){
    case Process::Type::data:
      datas_.push_back(move(hist));
      break;
    case Process::Type::background:
      backgrounds_.push_back(move(hist));
      break;
    case Process::Type::signal:
      signals_.push_back(move(hist));
      break;
    default:
      break;
    }
  }

  blank_.SetFillStyle(0);
  blank_.SetFillColor(kWhite);
  blank_.SetLineWidth(0);
  blank_.SetLineColor(kWhite);
  blank_.SetMarkerSize(0);
  blank_.SetMarkerColor(kWhite);
}

/*! \brief Produce and save formatted plots at given luminosity

  \param[in] luminosity The integrated luminosity with which to draw the plot
*/
void HistoStack::Print(double luminosity,
                       const string &subdir){
  luminosity_ = luminosity;
  for(const auto &opt: plot_options_){
    this_opt_ = opt;
    this_opt_.MakeSane();
    gStyle->SetColorModelPS(this_opt_.UseCMYK());
    gROOT->ForceStyle();
    RefreshScaledHistos();
    SetRanges();
    ApplyStyles();
    AdjustFillStyles();

    double bot_min, bot_max;
    vector<TH1D> bot_plots = GetBottomPlots(bot_min, bot_max);
    //I don't know why I can't make this in GetBottomPlots...
    TGraphAsymmErrors bottom_background;
    if(this_opt_.Bottom() != BottomType::off){
      bottom_background = TGraphAsymmErrors(&bot_plots.back());
      bottom_background.SetMinimum(this_opt_.RatioMinimum());
      bottom_background.SetMaximum(this_opt_.RatioMaximum());
      bot_plots.pop_back();
    }

    TGraphAsymmErrors bkg_error = GetBackgroundError();

    StripTopPlotLabels();
    TLine horizontal = GetBottomHorizontal();
    vector<TLine> cut_vals = GetCutLines(GetMinDraw(), GetMaxDraw(), true);
    vector<TLine> bot_cuts = GetCutLines(bot_min, bot_max, false);

    unique_ptr<TCanvas> full;
    unique_ptr<TPad> top, bottom;
    GetPads(full, top, bottom);

    if(this_opt_.AutoYAxis()) FixYAxis(bot_plots);

    if(this_opt_.Bottom() != BottomType::off){
      bottom->cd();

      string draw_opt = "e0";
      for(auto &h: bot_plots){
        h.Draw(draw_opt.c_str());
        draw_opt = "e0 same";
      }
      bottom_background.Draw("2 same");

      horizontal.Draw("same");

      bottom->RedrawAxis();
      bottom->RedrawAxis("g");
      for(auto &cut: bot_cuts) cut.Draw();
    }

    top->cd();
    if(this_opt_.YAxis() == YAxisType::log) top->SetLogy(true);
    else top->SetLogy(false);

    string draw_opt = "hist";
    DrawAll(backgrounds_, draw_opt);
    if(this_opt_.ShowBackgroundError() && backgrounds_.size()) bkg_error.Draw("2 same");
    DrawAll(signals_, draw_opt, true);
    ReplaceAll(draw_opt, "hist", "ep");
    DrawAll(datas_, draw_opt, true);
    for(auto &cut: cut_vals) cut.Draw();

    vector<shared_ptr<TLegend> > legends = GetLegends();
    for(auto &legend: legends){
      legend->Draw();
    }

    top->RedrawAxis();
    top->RedrawAxis("g");

    vector<shared_ptr<TLatex> > title_text = GetTitleTexts();
    for(auto &x: title_text){
      x->Draw();
    }

    // Printing values to terminal
    if(this_opt_.PrintVals()){
      TH1D *hdata = (datas_.size() ? &(datas_[0]->scaled_hist_) : 0);
      TH1D *hmc = (backgrounds_.size() ? &(backgrounds_[0]->scaled_hist_) : 0);
      TH1D *hbot = (bot_plots.size() ? &(bot_plots[0]) : 0);
      if(hdata==0 || hmc==0 || hbot==0) cout<<"Printing values not supported yet without Data, MC, or ratio"<<endl;
      else {
	int digits = floor(log10(max(hdata->GetBinContent(hdata->GetMaximumBin()), 
				     hmc->GetBinContent(hmc->GetMaximumBin())))+1.);
	//// Digits for error are calculated with the sqrt, and added 2 to print with one decimal
	int edigits = floor(log10(sqrt(max(hdata->GetMaximum(), hmc->GetMaximum())))+1.)+2;
	cout<<endl<<"Printing values for "<<definition_.Name()<<". Data/MC = "
	    <<RoundNumber(hdata->Integral(), 2,hmc->Integral()) <<endl;
	for(int bin=1; bin<=hdata->GetNbinsX(); bin++){
	  cout<<"Bin "<<setw(5)<<hdata->GetBinLowEdge(bin)<<","<<setw(5)<<hdata->GetBinLowEdge(bin+1)<<": Data = ";
	  cout<<setw(digits)<<hdata->GetBinContent(bin)<<" +- "<<setw(edigits)<<RoundNumber(hdata->GetBinError(bin),1);
	  cout<<", MC = "<<setw(digits+2)<<RoundNumber(hmc->GetBinContent(bin),1)<<" +- "
	      <<setw(edigits)<<RoundNumber(hmc->GetBinError(bin),1);
	  if(this_opt_.Bottom() != BottomType::off)
	    cout<<"   ->    Ratio  = "<<setw(6)<<RoundNumber(hbot->GetBinContent(bin),3)
		<<" +- "<<setw(5)<<RoundNumber(hbot->GetBinError(bin),3);
	  cout<<endl;
	} // Loop over histogram bins
      }
    }

    if(subdir != "") mkdir(("plots/"+subdir).c_str(), 0777);
    string base_name = subdir != ""
      ? "plots/"+subdir+"/"+definition_.Name()
      : "plots/"+definition_.Name();
    for(const auto &ext: this_opt_.FileExtensions()){
      string full_name = base_name+"__"+this_opt_.TypeString()+'.'+ext;
      full->Print(full_name.c_str());
      cout << "open " << full_name << endl;
    }
  }
}

set<const Process*> HistoStack::GetProcesses() const{
  set<const Process*> processes;
  for(const auto &proc: backgrounds_){
    processes.insert(proc->process_.get());
  }
  for(const auto &proc: signals_){
    processes.insert(proc->process_.get());
  }
  for(const auto &proc: datas_){
    processes.insert(proc->process_.get());
  }
  return processes;
}

Figure::FigureComponent * HistoStack::GetComponent(const Process *process){
  const auto &component_list = GetComponentList(process);
  for(const auto &component: component_list){
    if(component->process_.get() == process){
      return component.get();
    }
  }
  DBG("Could not find histogram for process "+process->name_+".");
  return nullptr;
}

/*!\brief Generates stacked and scaled histograms from unstacked and unscaled
  ones

  Sets bin contents for all required HistoStack::SingleHist::scaled_hist_ to the
  appropriate values using the HistoStack::SingleHist::raw_hist_ containing the
  unstacked contents at 1 fb^{-1}
*/
void HistoStack::RefreshScaledHistos(){
  InitializeHistos();
  MergeOverflow();
  ScaleHistos();
  StackHistos();
  NormalizeHistos();
}

/*!\brief Sets all HistoStack::SingleHist::scaled_hist_ to corresponding
  HistoStack::SingleHist::raw_hist_
*/
void HistoStack::InitializeHistos() const{
  for(auto &hist: backgrounds_){
    hist->scaled_hist_ = hist->raw_hist_;
    hist->scaled_hist_.SetName(("bkg_"+hist->process_->name_+"_"+counter()).c_str());
  }
  for(auto &hist: signals_){
    hist->scaled_hist_ = hist->raw_hist_;
    hist->scaled_hist_.SetName(("sig_"+hist->process_->name_+"_"+counter()).c_str());
  }
  for(auto &hist: datas_){
    hist->scaled_hist_ = hist->raw_hist_;
    hist->scaled_hist_.SetName(("dat_"+hist->process_->name_+"_"+counter()).c_str());
  }
}

/*!\brief Moves overflow (underflow) contents into last (first) visible bin
 */
void HistoStack::MergeOverflow() const{
  bool underflow = false, overflow = false;
  switch(this_opt_.Overflow()){
  default:
  case OverflowType::none:
    underflow = false;
    overflow = false;
    break;
  case OverflowType::underflow:
    underflow = true;
    overflow = false;
    break;
  case OverflowType::overflow:
    underflow = false;
    overflow = true;
    break;
  case OverflowType::both:
    underflow = true;
    overflow = true;
    break;
  }

  for(auto &hist: backgrounds_){
    ::MergeOverflow(hist->scaled_hist_, underflow, overflow);
  }
  for(auto &hist: signals_){
    ::MergeOverflow(hist->scaled_hist_, underflow, overflow);
  }
  for(auto &hist: datas_){
    ::MergeOverflow(hist->scaled_hist_, underflow, overflow);
  }
}

/*!\brief Scales histograms to required luminosity
 */
void HistoStack::ScaleHistos() const{
  for(auto &hist: backgrounds_){
    AdjustDensityForBinWidth(hist->scaled_hist_);
    hist->scaled_hist_.Scale(luminosity_);
  }
  for(auto &hist: signals_){
    AdjustDensityForBinWidth(hist->scaled_hist_);
    hist->scaled_hist_.Scale(luminosity_);
  }
  for(auto &hist: datas_){
    AdjustDensityForBinWidth(hist->scaled_hist_);
  }
}

/*!\brief Stacks histograms if necessary for current plot style
 */
void HistoStack::StackHistos() const{
  if(this_opt_.Stack() == StackType::signal_overlay
     || this_opt_.Stack() == StackType::signal_on_top
     || this_opt_.Stack() == StackType::data_norm){
    for(size_t ibkg = backgrounds_.size() - 2; ibkg < backgrounds_.size(); --ibkg){
      backgrounds_.at(ibkg)->scaled_hist_ = backgrounds_.at(ibkg)->scaled_hist_ + backgrounds_.at(ibkg+1)->scaled_hist_;
    }
    if(backgrounds_.size() && this_opt_.Stack() == StackType::signal_on_top){
      for(auto &hist: signals_){
        hist->scaled_hist_ = hist->scaled_hist_ + backgrounds_.front()->scaled_hist_;
      }
    }
  }
}

/*!\brief Normalize histograms to data or 100%*(bin width) if needed for current
  style
*/
void HistoStack::NormalizeHistos() const{
  mc_scale_ = 1.;
  mc_scale_error_ = 1.;
  if(this_opt_.Stack() == StackType::data_norm){
    if(datas_.size() == 0 || backgrounds_.size() == 0) return;
    int nbins = definition_.Nbins();
    double data_error, mc_error;
    double data_norm = datas_.front()->scaled_hist_.IntegralAndError(0, nbins+1, data_error, "width");
    double mc_norm = backgrounds_.front()->scaled_hist_.IntegralAndError(0, nbins+1, mc_error, "width");
    mc_scale_ = data_norm/mc_norm;
    mc_scale_error_ = hypot(data_norm*mc_error, mc_norm*data_error)/(mc_norm*mc_norm);
    for(auto &hist: backgrounds_){
      hist->scaled_hist_.Scale(mc_scale_);
    }
    for(auto h = datas_.begin(); h != datas_.end(); ++h){
      auto &hist = *h;
      double dumb;
      double this_integral = hist->scaled_hist_.IntegralAndError(0, nbins+1, dumb, "width");
      hist->scaled_hist_.Scale(this_integral == 0. ? 1. : data_norm/this_integral);
    }
  }else if(this_opt_.Stack() == StackType::shapes){
    for(auto &hist: backgrounds_){
      Normalize(hist->scaled_hist_, 100., true);
    }
    for(auto &hist: signals_){
      Normalize(hist->scaled_hist_, 100., true);
    }
    for(auto &hist: datas_){
      Normalize(hist->scaled_hist_, 100., true);
    }
  }
}

/*!\brief Set y-axis plotting range
 */
void HistoStack::SetRanges() const{
  double the_min = GetMinDraw();
  double the_max = GetMaxDraw();

  double ratio = GetLegendRatio();

  double top, bottom;
  switch(this_opt_.YAxis()){
  default:
  case YAxisType::linear:
    bottom = the_min >= 0. ? 0. : the_min;
    top = bottom+ratio*(the_max-bottom);
    break;
  case YAxisType::log:
    bottom = the_min > this_opt_.LogMinimum() ? the_min : this_opt_.LogMinimum();
    top = exp(log(bottom)+ratio*(log(the_max)-log(bottom)));
    break;
  }

  for(auto &hist: backgrounds_){
    hist->scaled_hist_.SetMinimum(bottom);
    hist->scaled_hist_.SetMaximum(top);
  }
  for(auto &hist: signals_){
    hist->scaled_hist_.SetMinimum(bottom);
    hist->scaled_hist_.SetMaximum(top);
  }
  for(auto &hist: datas_){
    hist->scaled_hist_.SetMinimum(bottom);
    hist->scaled_hist_.SetMaximum(top);
  }
}

/*!\brief Set label styles and title for all histograms
 */
void HistoStack::ApplyStyles() const{
  for(auto &hist: backgrounds_){
    StyleHisto(hist->scaled_hist_);
  }
  for(auto &hist: signals_){
    StyleHisto(hist->scaled_hist_);
  }
  for(auto &hist: datas_){
    StyleHisto(hist->scaled_hist_);
  }
}

/*!\brief Set label styles and title for a histogram

  \param[in,out] h Histogram to be adjusted
*/
void HistoStack::StyleHisto(TH1D &h) const{
  h.GetXaxis()->SetTitleOffset(this_opt_.XTitleOffset());
  h.GetYaxis()->SetTitleOffset(this_opt_.YTitleOffset());
  h.SetNdivisions(this_opt_.NDivisions(), "xyz");
  h.SetLabelSize(this_opt_.LabelSize(), "xyz");
  h.SetTitleSize(this_opt_.TitleSize(), "xyz");
  h.SetLabelFont(this_opt_.Font(), "xyz");
  h.SetTitleFont(this_opt_.Font(), "xyz");
  
  double bin_width = (definition_.Bins().back()-definition_.Bins().front())/(definition_.Nbins());
  
  ostringstream title;
  switch(this_opt_.Stack()){
  default:
    DBG("Unrecognized stack option " << static_cast<int>(this_opt_.Stack()) << ".");
  case StackType::signal_overlay:
  case StackType::signal_on_top:
  case StackType::data_norm:
  case StackType::lumi_shapes:
    if(definition_.units_ == "" && bin_width == 1){
      title << "Entries";    
      break;
    }
    else{
      title << "Entries/(" << bin_width;
      if(definition_.units_ != "") title << " " << definition_.units_;
      title << ")";
      break;
    }
  case StackType::shapes:
    if(definition_.units_ == "" && bin_width == 1){
      title << "% entries";
      break;
    }
    else{
      title << "% entries/(" << bin_width;
      if(definition_.units_ != "") title << " " << definition_.units_;
      title << ")";
      break;
    }
  }
  
  h.GetYaxis()->SetTitle(title.str().c_str());
}

/*!\brief Make histograms a hollow line for unstacked styles
 */
void HistoStack::AdjustFillStyles() const{
  if(this_opt_.BackgroundsStacked()) return;

  for(auto &bkg: backgrounds_){
    TH1D &h = bkg->scaled_hist_;
    h.SetFillStyle(0);
    h.SetLineColor(h.GetFillColor());
    h.SetLineWidth(5);
  }
}

/*!\brief Generated canvas and pads for top and bottom plots

  To make object positioning simpler, the top and bottom pads span the entire
  canvas and only the margins are adjusted. This has the added bonus of making
  the "0" on the top y-axis visible for free.

  \param[out] c Full plot canvas

  \param[out] top Pad for main plot

  \param[out] bottom Pad for bottom plot (e.g. ratio plot)
*/
void HistoStack::GetPads(unique_ptr<TCanvas> &c,
                         unique_ptr<TPad> &top,
                         unique_ptr<TPad> &bottom) const{
  c.reset(new TCanvas(("canvas_"+counter()).c_str(), "canvas", this_opt_.CanvasWidth(),
                      this_opt_.CanvasHeight()));
  c->cd();
  top.reset(new TPad(("top_pad_"+counter()).c_str(), "top_pad", 0., 0., 1., 1.));
  bottom.reset(new TPad(("bottom_pad_"+counter()).c_str(), "bottom_pad", 0., 0., 1., 1.));
  c->SetMargin(0., 0., 0., 0.);
  c->SetTicks(1,1);
  c->SetFillStyle(4000);
  top->SetTicks(1,1);
  top->SetFillStyle(4000);
  bottom->SetTicks(1,1);
  bottom->SetFillStyle(4000);
  if(this_opt_.Bottom() == BottomType::off){
    top->SetMargin(this_opt_.LeftMargin(),
                   this_opt_.RightMargin(),
                   this_opt_.BottomMargin(),
                   this_opt_.TopMargin());
  }else{
    top->SetMargin(this_opt_.LeftMargin(),
                   this_opt_.RightMargin(),
                   this_opt_.BottomHeight(),
                   this_opt_.TopMargin());
    bottom->SetMargin(this_opt_.LeftMargin(),
                      this_opt_.RightMargin(),
                      this_opt_.BottomMargin(),
                      1.-this_opt_.BottomHeight());
  }
  bottom->Draw();
  top->Draw();
}

/*!\brief Adjust y-axis title offset based on y-axis range

  \param[in,out] bottom_plots Ratio or other bottom-half plots whose title
  offsets also need to be adjusted
*/
void HistoStack::FixYAxis(vector<TH1D> &bottom_plots) const{
  double offset = this_opt_.YTitleOffset();
  if(this_opt_.YAxis() == YAxisType::log){
    offset = 1.5;
  }else{
    double the_max = GetMaxDraw()*GetLegendRatio();
    int digits = fabs(floor(log10(the_max))-1)+2;
    digits = max(2, min(6, digits));//For huge axis scale, reverts to scientific notation

    //Scale offset by good empirical numbers
    offset = 0.6+0.25*digits;
  
}
  for(auto &hist: backgrounds_){
    hist->scaled_hist_.SetTitleOffset(offset, "y");
  }
  for(auto &hist: signals_){
    hist->scaled_hist_.SetTitleOffset(offset, "y");
  }
  for(auto &hist: datas_){
    hist->scaled_hist_.SetTitleOffset(offset, "y");
  }
  for(auto &hist: bottom_plots){
    if(datas_.size()>0)hist.SetTitleOffset(offset, "y");
  }
}

/*!\brief Get text to print at top of plot

  Depending on current plot style, this may be the CMS {Preliminary, Simulation,
  etc...} with luminosity or the cut and weight for the plot

  \return List of text items to be printed in title region of plot
*/
vector<shared_ptr<TLatex> > HistoStack::GetTitleTexts() const{
  vector<shared_ptr<TLatex> > out;
  double left = this_opt_.LeftMargin();
  double right = 1.-this_opt_.RightMargin();
  double bottom = 1.-this_opt_.TopMargin();
  double top = 1.;
  if(this_opt_.Title() == TitleType::info){
    if(definition_.Title() != ""){
      out.push_back(make_shared<TLatex>(0.5*(left+right), 0.5*(bottom+top),
                                        definition_.Title().c_str()));
      out.back()->SetNDC();
      out.back()->SetTextAlign(22);
      out.back()->SetTextFont(this_opt_.Font());

      //Adjust title to fit in available space
      double max_width, max_height;
      GetTitleSize(max_width, max_height, true);
      UInt_t width, height;
      out.back()->GetBoundingBox(width, height);
      while(width > max_width || height > max_height){
        out.back()->SetTextSize(0.8*out.back()->GetTextSize());
        out.back()->GetBoundingBox(width, height);
      }
      while(width < 0.5*max_width && height < 0.5*max_height){
        out.back()->SetTextSize(1.25*out.back()->GetTextSize());
        out.back()->GetBoundingBox(width, height);
      }
    }
  }else{
    string extra;
    switch(this_opt_.Title()){
    case TitleType::preliminary: extra = "Preliminary"; break;
    case TitleType::simulation: extra = "Simulation"; break;
    case TitleType::supplementary: extra = "Supplementary"; break;
    case TitleType::data: extra = ""; break;
    case TitleType::info:
    default:
      ERROR("Did not understand title type "+to_string(static_cast<int>(this_opt_.Title())));
    }
    out.push_back(make_shared<TLatex>(left, 0.5*(bottom+top),
                                      ("#font[62]{CMS}#scale[0.76]{#font[52]{ "+extra+"}}").c_str()));
    out.back()->SetNDC();
    out.back()->SetTextAlign(12);
    out.back()->SetTextFont(this_opt_.Font());

    ostringstream oss;
    oss << luminosity_ << " fb^{-1} (13 TeV)" << flush;
    out.push_back(make_shared<TLatex>(right, 0.5*(bottom+top),
                                      oss.str().c_str()));
    out.back()->SetNDC();
    out.back()->SetTextAlign(32);
    out.back()->SetTextFont(this_opt_.Font());
  }
  return out;
}

/*!\brief Get uncertainty on total background

  \return Graph representing the MC background error band
*/
TGraphAsymmErrors HistoStack::GetBackgroundError() const{
  TGraphAsymmErrors g;
  if(backgrounds_.size() == 0){
    TH1D h("", "", definition_.Nbins(), &definition_.Bins().at(0));
    g = TGraphAsymmErrors(&h);
  }else{
    g = TGraphAsymmErrors(&(backgrounds_.front()->scaled_hist_));
  }
  g.SetFillStyle(3002);
  // set the color of the error band to the line color, accomodating data-to-data plots
  g.SetFillColor(backgrounds_.front()->scaled_hist_.GetLineColor());
  g.SetLineWidth(0);
  g.SetMarkerSize(0);
  return g;
}

/*!\brief Get vertical lines at cut values

  \param[in] y_min Lower bound of y-axis

  \param[in] y_max Upper bound of y-axis

  \return Lines at x-coordinate of cut value and y-coordinates running from
  bottom of plot to bottom of legend
*/
vector<TLine> HistoStack::GetCutLines(double y_min, double y_max, bool adjust_bottom) const{
  double bottom = y_min;
  if(adjust_bottom){
    switch(this_opt_.YAxis()){
    default:
      DBG("Bad YAxis type " << static_cast<int>(this_opt_.YAxis()));
    case YAxisType::linear: bottom = y_min >= 0. ? 0. : y_min; break;
    case YAxisType::log:    bottom = y_min > this_opt_.LogMinimum() ? y_min : this_opt_.LogMinimum(); break;
    }
  }
  vector<TLine> out(definition_.cut_vals_.size());
  for(double cut: definition_.cut_vals_){
    out.emplace_back(cut, bottom, cut, y_max);
    out.back().SetNDC(false);
    out.back().SetLineStyle(2);
    out.back().SetLineColor(kBlack);
    out.back().SetLineWidth(3);
  }

  return out;
}

/*!\brief Get ratio or other plots drawn on the lower pad

  \param [out] the_min Y-axis minimum across plots for lower pad

  \param [out] the_max Y-axis maximum across plots for lower pad

  \return Set of plots to be drawn on lower pad. These may be ratio plots or
  something else depending on the current plot style
*/
std::vector<TH1D> HistoStack::GetBottomPlots(double &the_min, double &the_max) const{
  if(this_opt_.Bottom() == BottomType::off) return vector<TH1D>();

  TH1D denom;
  vector<TH1D> out;
  if(backgrounds_.size() != 0){
    denom = backgrounds_.front()->scaled_hist_;
  }else if(datas_.size() != 0){
    denom = datas_.front()->scaled_hist_;
  }else if(signals_.size() != 0){
    denom = signals_.front()->scaled_hist_;
  }else{
    ERROR("No histograms available to make bottom plot");
  }
  bool stacked;
  switch(this_opt_.Stack()){
  case StackType::signal_overlay:
  case StackType::signal_on_top:
  case StackType::data_norm:
    stacked = true; break;
  case StackType::lumi_shapes:
  case StackType::shapes:
    stacked = false; break;
  default:
    ERROR("Bad stack type: "+to_string(static_cast<int>(this_opt_.Stack())));
    break;
  }
  if(stacked && backgrounds_.size()){
    out.push_back(backgrounds_.front()->scaled_hist_);
    out.back().SetName(("bot_plot_bkg_"+backgrounds_.front()->process_->name_+"_"+counter()).c_str());
  }else{
    for(const auto &h: backgrounds_){
      out.push_back(h->scaled_hist_);
      out.back().SetName(("bot_plot_bkg_"+h->process_->name_+"_"+counter()).c_str());
    }
  }
  for(const auto &h: datas_){
    out.push_back(h->scaled_hist_);
    out.back().SetName(("bot_plot_data_"+h->process_->name_+"_"+counter()).c_str());
  }
  if(!stacked){
    for(const auto &h: signals_){
      out.push_back(h->scaled_hist_);
      out.back().SetName(("bot_plot_sig_"+h->process_->name_+"_"+counter()).c_str());
    }
  }
  if(!out.size()) return vector<TH1D>();
  TH1D band = out.front();
  for(size_t i = 0; (i+1) < out.size(); ++i){
    out.at(i) = out.at(i+1);
  }
  out.back() = band;
  out.back().SetFillStyle(3002);
  // makes ratio error band colored for data-to-data plots as well
  if (backgrounds_.front()->scaled_hist_.GetLineColor()!=kBlack)
    out.back().SetFillColor(backgrounds_.front()->scaled_hist_.GetLineColor());
  out.back().SetLineWidth(0);
  out.back().SetMarkerStyle(0);
  out.back().SetMarkerSize(0);
  out.back().SetName(("bot_plot_band_"+counter()).c_str());

  for(int bin = 0; bin <= denom.GetNbinsX()+1; ++bin){
    denom.SetBinError(bin, 0.);
  }

  switch(this_opt_.Bottom()){
  case BottomType::ratio:
    for(auto &h: out){
      h.Divide(&denom);
    }
    break;
  case BottomType::diff:
    for(auto &h: out){
      h = h - denom;
    }
    break;
  case BottomType::off:
  default:
    ERROR("Bad type for bottom plot: "+to_string(static_cast<int>(this_opt_.Bottom())));
    break;
  }

  the_min = numeric_limits<double>::infinity();
  the_max = -numeric_limits<double>::infinity();
  for(auto &h: out){
    h.SetNdivisions(this_opt_.NDivisionsBottom(), "y");
    for(int bin = 1; bin <= h.GetNbinsX(); ++bin){
      double hi = h.GetBinContent(bin)+h.GetBinErrorUp(bin);
      double lo = h.GetBinContent(bin)-fabs(h.GetBinErrorLow(bin));
      if(hi>the_max) the_max = hi;
      if(lo<the_min) the_min = lo;
    }
  }

  if(this_opt_.Bottom() == BottomType::ratio){
    the_min = this_opt_.RatioMinimum();
    the_max = this_opt_.RatioMaximum();
    for(auto &h: out){
      if(datas_.size() != 0) h.GetYaxis()->SetTitle("Data/MC");
      else {
        string label = "#frac{MC}{"+backgrounds_.front()->process_->name_+"}";
        h.GetYaxis()->SetTitle(label.c_str());
        h.SetTitleSize(h.GetTitleSize("y")/1.5,"y");
        h.SetTitleOffset(2.1,"y");
      }
      h.SetMinimum(the_min);
      h.SetMaximum(the_max);
    }
  }else if(this_opt_.Bottom() == BottomType::diff){
    for(auto &h: out){
      if(datas_.size() != 0) h.GetYaxis()->SetTitle("Data-MC");
      //else h.GetYaxis()->SetTitle(backgrounds_.front()->process_->name_.c_str());
      h.SetMinimum(the_min);
      h.SetMaximum(the_max);
    }
}
  return out;
}

/*!\brief Get horizontal line drawn at "agreement" value for bottom plots

  E.g. Line is at 1 for ratio plots, 0 for difference plots, etc.

  \return Line at appropriate height depending on plot style
*/
TLine HistoStack::GetBottomHorizontal() const{
  double left = definition_.Bins().front();
  double right = definition_.Bins().back();
  double y;
  switch(this_opt_.Bottom()){
  case BottomType::ratio: y = 1.; break;
  case BottomType::diff: y = 0.; break;
  case BottomType::off: y = 0.; break;
  default:
    y = 0.;
    DBG("Invalid BottomType: " << to_string(static_cast<int>(this_opt_.Bottom())));
  }

  TLine line(left, y, right, y);
  line.SetNDC(false);
  line.SetLineStyle(2);
  line.SetLineColor(kBlack);
  line.SetLineWidth(2);
  return line;
}

/*!Remove x-axis labels and title from plots in top pad if necessary
 */
void HistoStack::StripTopPlotLabels() const{
  if(this_opt_.Bottom() == BottomType::off) return;
  for(auto &hist: backgrounds_){
    StripXLabels(hist->scaled_hist_);
  }
  for(auto &hist: signals_){
    StripXLabels(hist->scaled_hist_);
  }
  for(auto &hist: datas_){
    StripXLabels(hist->scaled_hist_);
  }
}

/*!\brief Get highest drawn point below max_bound across all component
  histograms

  \param[in] max_bound Only consider points below this value in finding the
  maximum

  \return The highest drawn point below max_bound across all component
  histograms
*/
double HistoStack::GetMaxDraw(double max_bound) const{
  double the_max = -numeric_limits<double>::infinity();
  for(const auto &hist: backgrounds_){
    double this_max = hist->GetMax(max_bound, this_opt_.ShowBackgroundError());
    if(this_max > the_max && this_max < max_bound){
      the_max = this_max;
    }
  }
  for(const auto &hist: signals_){
    double this_max = hist->GetMax(max_bound, false);
    if(this_max > the_max && this_max < max_bound){
      the_max = this_max;
    }
  }
  for(const auto &hist: datas_){
    double this_max = hist->GetMax(max_bound, true);
    if(this_max > the_max && this_max < max_bound){
      the_max = this_max;
    }
  }
  return the_max;
}

/*!\brief Get lowest drawn point above min_bound across all component histograms

  \param[in] min_bound Only consider points above this value in finding the
  minimum

  \return The lowest drawn point greater than min_bound across all component
  histograms
*/
double HistoStack::GetMinDraw(double min_bound) const{
  double the_min = numeric_limits<double>::infinity();
  for(const auto &hist: backgrounds_){
    double this_min = hist->GetMin(min_bound, this_opt_.ShowBackgroundError());
    if(this_min < the_min && this_min > min_bound){
      the_min = this_min;
    }
  }
  for(const auto &hist: signals_){
    double this_min = hist->GetMin(min_bound, false);
    if(this_min < the_min && this_min > min_bound){
      the_min = this_min;
    }
  }
  for(const auto &hist: datas_){
    double this_min = hist->GetMin(min_bound, true);
    if(this_min < the_min && this_min > min_bound){
      the_min = this_min;
    }
  }
  return the_min;
}

/*!\brief Get list of legends emulating single legend with multiple columns

  Legends are filled down-columns first, then across rows. Data samples are
  added first, then signals, then backgrounds, with the order within each group
  preserved from HistoStack::HistoStack()

  \return Legends with all processes and possibly MC normalization if plot style
  requires it
*/
vector<shared_ptr<TLegend> > HistoStack::GetLegends(){
  size_t n_entries = datas_.size() + signals_.size() + backgrounds_.size();
  if(this_opt_.DisplayLumiEntry()) ++n_entries;
  size_t n_columns = min(n_entries, static_cast<size_t>(this_opt_.LegendColumns()));

  double left = this_opt_.LeftMargin()+this_opt_.LegendPad();
  double top = 1.-this_opt_.TopMargin()-this_opt_.LegendPad();
  double bottom = top-this_opt_.TrueLegendHeight(n_entries);

  double delta_x = this_opt_.TrueLegendWidth(n_entries);
  vector<shared_ptr<TLegend> > legends(n_columns);
  for(size_t i = 0; i < n_columns; ++i){
    double x = left+i*delta_x;
    legends.at(i) = make_shared<TLegend>(x, bottom, x+this_opt_.LegendMarkerWidth(), top);
    legends.at(i)->SetFillStyle(0);
    legends.at(i)->SetBorderSize(0);
    legends.at(i)->SetTextSize(this_opt_.TrueLegendEntryHeight(n_entries));
    legends.at(i)->SetTextFont(this_opt_.Font());
  }

  size_t entries_added = 0;
  AddEntries(legends, datas_, "lep", n_entries, entries_added);
  AddEntries(legends, signals_, "l", n_entries, entries_added);
  AddEntries(legends, backgrounds_, this_opt_.BackgroundsStacked() ? "f" : "l", n_entries, entries_added);

  //Add a dummy legend entry to display MC normalization
  if(this_opt_.DisplayLumiEntry()){
    auto &leg = legends.at(GetLegendIndex(entries_added, n_entries, legends.size()));
    ostringstream label;
    label << fixed << setprecision(1) << "L=" << luminosity_ << " fb^{-1}";
    if(this_opt_.Stack() == StackType::data_norm && datas_.size() > 0){
      label << ", (" << 100.*mc_scale_ << "#pm" << 100.*mc_scale_error_ << ")%";
    }
    auto entry = leg->AddEntry(&blank_, label.str().c_str(), "f");
    entry->SetFillStyle(0);
    entry->SetFillColor(kWhite);
    entry->SetLineWidth(0);
    entry->SetLineColor(kWhite);
    entry->SetMarkerStyle(0);
    entry->SetMarkerStyle(kWhite);
  }

  return legends;
}

/*!\brief Distribute processes from list of histograms across legends

  \param[in,out] legends Legends to which entries are added

  \param[in] hists Component histograms whose processes go in legends

  \param[in] style Option string specifying legend marker type. Subset of "flep"

  \param[in] n_entries Number of entries that need to fit in legends

  \param[in,out] entries_added Entries already in legend
*/
void HistoStack::AddEntries(vector<shared_ptr<TLegend> > &legends,
                            const vector<unique_ptr<SingleHist> > &hists,
                            const string &style,
                            size_t n_entries,
                            size_t &entries_added) const{
  for(auto h = hists.cbegin(); h != hists.cend(); ++h){
    size_t legend_index = GetLegendIndex(entries_added, n_entries, legends.size());
    TLegend &legend = *legends.at(legend_index);
    string label = (*h)->process_->name_.c_str();
    if(this_opt_.Title() == TitleType::info){
      double value;
      switch(this_opt_.Stack()){
      default:
        DBG("Bad stack option: " << static_cast<int>(this_opt_.Stack()));
      case StackType::signal_overlay:
      case StackType::signal_on_top:
      case StackType::data_norm:
        value = GetYield(h);
        if(value>=1.){
          label += " [N=" + FixedDigits(value, 2) + "]";
        }else{
          label += " [N=" + FixedDigits(value, 1) + "]";
        }
        break;
      case StackType::lumi_shapes:
      case StackType::shapes:
        value = GetMean(h);
        label += " [#mu=" + FixedDigits(value,3) + "]";
        break;
      }
    }
    //Shrink text size if label is long
    double fudge_factor = 0.25;//Not sure how TLegend width affects marker width, but this seems to work
    double max_width = (this_opt_.TrueLegendWidth(n_entries)-this_opt_.LegendMarkerWidth()*fudge_factor) * this_opt_.CanvasWidth();
    double max_height = this_opt_.TrueLegendEntryHeight(n_entries) * this_opt_.LegendDensity() * this_opt_.CanvasHeight();
    TLatex latex(0.5, 0.5, label.c_str());
    latex.SetTextSize(legend.GetTextSize());
    latex.SetTextFont(legend.GetTextFont());
    UInt_t width, height;
    latex.GetBoundingBox(width, height);
    while(width > max_width || height > max_height){
      latex.SetTextSize(0.95*latex.GetTextSize());
      for(auto &leg: legends){
        leg->SetTextSize(0.95*leg->GetTextSize());
      }
      latex.GetBoundingBox(width, height);
    }

    legend.AddEntry(&((*h)->scaled_hist_), label.c_str(), style.c_str());
    ++entries_added;
  }
}

/*!\brief Get factor by which to expand y-axis range to fit legend

  \return Factor by which the upper bound of the top plot's y-axis needs to be
  expanded to make room for the legend
*/
double HistoStack::GetLegendRatio() const{
  size_t num_plots = backgrounds_.size() + signals_.size() + datas_.size();
  if(this_opt_.DisplayLumiEntry()) ++num_plots;
  double legend_height = this_opt_.TrueLegendHeight(num_plots);
  double top_plot_height;
  if(this_opt_.Bottom() == BottomType::off){
    top_plot_height = 1.-this_opt_.TopMargin()-this_opt_.BottomMargin();
  }else{
    top_plot_height = 1.-this_opt_.TopMargin()-this_opt_.BottomHeight();
  }
  return top_plot_height/(top_plot_height-legend_height-2.*this_opt_.LegendPad());
}

/*!\brief Get integrated number of weighted entries in histogram

  Possibly varying bin widths are accounted for.

  \param[in] h Iterator to histogram in one of HistoStack::datas_,
  HistoStack::signals_, HistoStack::backgrounds_ for which total yield will be
  found
*/
double HistoStack::GetYield(std::vector<std::unique_ptr<SingleHist> >::const_iterator h) const{
  TH1D hist = (*h)->scaled_hist_;

  //Subtract underlying histogram
  if((*h)->process_->type_ == Process::Type::background
     && h != (--backgrounds_.cend())
     && this_opt_.BackgroundsStacked()){
    hist = hist - (*(++h))->scaled_hist_;
  }

  //Want yield, not area, so divide out average bin width
  //N.B.: Can't just use hist.Integral() in case of varying bin width
  double raw_integral = hist.Integral("width");
  int nbins = hist.GetNbinsX();
  double left = hist.GetBinLowEdge(1);
  double right = hist.GetBinLowEdge(nbins+1);
  double width = (right-left)/nbins;
  return raw_integral/width;
}

/*!\brief Get mean of histogram

  Possibly varying bin widths are accounted for.

  \param[in] h Iterator to histogram in one of HistoStack::datas_,
  HistoStack::signals_, HistoStack::backgrounds_ for which mean will be found
*/
double HistoStack::GetMean(std::vector<std::unique_ptr<SingleHist> >::const_iterator h) const{
  TH1D hist = (*h)->scaled_hist_;

  //Subtract underlying histogram
  if((*h)->process_->type_ == Process::Type::background
     && h != (--backgrounds_.cend())
     && this_opt_.BackgroundsStacked()){
    hist = hist - (*(++h))->scaled_hist_;
  }

  return hist.GetMean();
}

/*!\brief Get width and height of title region

  \param[out] width Width of title region

  \param[out] height Height of title region

  \param[in] in_pixels If true, measure dimensions in pixels. If false, measure
  in NDC.
*/
void HistoStack::GetTitleSize(double &width, double &height, bool in_pixels) const{
  width = 1.-this_opt_.LeftMargin()-this_opt_.RightMargin();
  height = this_opt_.TopMargin();
  if(in_pixels){
    width *= this_opt_.CanvasWidth();
    height *= this_opt_.CanvasHeight();
  }
}

const vector<unique_ptr<HistoStack::SingleHist> >& HistoStack::GetComponentList(const Process *process){
  switch(process->type_){
  case Process::Type::data:
    return datas_;
  case Process::Type::background:
    return backgrounds_;
  case Process::Type::signal:
    return signals_;
  default:
    ERROR("Did not understand process type "+to_string(static_cast<long>(process->type_))+".");
    return backgrounds_;
  }
}
