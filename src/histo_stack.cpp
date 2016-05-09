#include "histo_stack.hpp"

#include <cmath>

#include <algorithm>
#include <sstream>

#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMath.h"

#include "utilities.hpp"

using namespace std;
using namespace PlotOptTypes;

namespace{
  void DrawAll(vector<HistoStack::SingleHist> &hists,
               string &draw_opt){
    for(auto hist = hists.rbegin(); hist!= hists.rend(); ++hist){
      hist->scaled_hist_.Draw(draw_opt.c_str());
      draw_opt = "hist same";
    }
  }

  void StripXLabels(TH1D &h){
    h.GetXaxis()->SetTitle("");
    h.SetLabelSize(0., "x");
    h.SetTitleSize(0., "x");
  }

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

HistoStack::HistoStack(const vector<shared_ptr<Process> > &processes,
                       const HistoDef &definition,
                       const PlotOpt &plot_options):
  backgrounds_(),
  signals_(),
  datas_(),
  definition_(definition),
  plot_options_(plot_options){
  string x_title = definition.x_title_;
  if(definition.units_ != "") x_title += " ["+definition.units_+"]";

  double bin_width = (definition.Bins().back()-definition.Bins().front())/(definition.Nbins());

  ostringstream oss;
  switch(plot_options.Stack()){
  default:
    DBG("Unrecognized stack option " << static_cast<int>(plot_options.Stack()) << ".");
  case StackType::signal_overlay:
  case StackType::signal_on_top:
  case StackType::lumi_shapes:
    oss << "Entries/(" << bin_width;
    if(definition.units_ != "") oss << " " << definition.units_;
    oss << ")";
    break;
  case StackType::shapes:
    oss << "% entries/(" << bin_width;
    if(definition.units_ != "") oss << " " << definition.units_;
    oss << ")";
    break;
  }

  TH1D empty("", (";"+x_title+";"+oss.str()).c_str(), definition.Nbins(), &definition.Bins().at(0));
  empty.SetStats(false);
  empty.GetXaxis()->SetTitleOffset(plot_options_.XTitleOffset());
  empty.GetYaxis()->SetTitleOffset(plot_options_.YTitleOffset());
  empty.SetNdivisions(plot_options_.NDivisions(), "xyz");
  empty.SetLabelSize(plot_options_.LabelSize(), "xyz");
  empty.SetTitleSize(plot_options_.TitleSize(), "xyz");
  empty.SetLabelFont(plot_options_.Font(), "xyz");
  empty.SetTitleFont(plot_options_.Font(), "xyz");
  empty.Sumw2(true);
  for(const auto &process: processes){
    SingleHist hist(process, empty);
    hist.raw_hist_.SetFillColor(process->GetFillColor());
    hist.raw_hist_.SetFillStyle(process->GetFillStyle());
    hist.raw_hist_.SetLineColor(process->GetLineColor());
    hist.raw_hist_.SetLineStyle(process->GetLineStyle());
    hist.raw_hist_.SetLineWidth(process->GetLineWidth());
    hist.raw_hist_.SetMarkerColor(process->GetMarkerColor());
    hist.raw_hist_.SetMarkerStyle(process->GetMarkerStyle());
    hist.raw_hist_.SetMarkerSize(process->GetMarkerSize());

    switch(process->type_){
    case Process::Type::data:
      datas_.push_back(hist);
      break;
    case Process::Type::background:
      backgrounds_.push_back(hist);
      break;
    case Process::Type::signal:
      signals_.push_back(hist);
      break;
    default:
      break;
    }
  }
  }

void HistoStack::StripTopPlotLabels(){
  if(plot_options_.Bottom() == BottomType::off) return;
  for(auto &hist: backgrounds_){
    StripXLabels(hist.scaled_hist_);
  }
  for(auto &hist: signals_){
    StripXLabels(hist.scaled_hist_);
  }
  for(auto &hist: datas_){
    StripXLabels(hist.scaled_hist_);
  }
}

void HistoStack::GetPads(unique_ptr<TCanvas> &c,
                         unique_ptr<TPad> &top,
                         unique_ptr<TPad> &bottom) const{
  c.reset(new TCanvas("", "", plot_options_.CanvasWidth(),
                      plot_options_.CanvasHeight()));
  c->cd();
  top.reset(new TPad("", "", 0., 0., 1., 1.));
  bottom.reset(new TPad("", "", 0., 0., 1., 1.));
  c->SetMargin(0., 0., 0., 0.);
  c->SetTicks(1,1);
  c->SetFillStyle(4000);
  top->SetTicks(1,1);
  top->SetFillStyle(4000);
  bottom->SetTicks(1,1);
  bottom->SetFillStyle(4000);
  if(plot_options_.Bottom() == BottomType::off){
    top->SetMargin(plot_options_.LeftMargin(),
                   plot_options_.RightMargin(),
                   plot_options_.BottomMargin(),
                   plot_options_.TopMargin());
  }else{
    top->SetMargin(plot_options_.LeftMargin(),
                   plot_options_.RightMargin(),
                   plot_options_.BottomHeight(),
                   plot_options_.TopMargin());
    bottom->SetMargin(plot_options_.LeftMargin(),
                      plot_options_.RightMargin(),
                      plot_options_.BottomMargin(),
                      1.-plot_options_.BottomHeight());
  }
  bottom->Draw();
  top->Draw();
}

vector<shared_ptr<TLatex> > HistoStack::GetTitleTexts(double luminosity, double left_margin) const{
  vector<shared_ptr<TLatex> > out;
  double left = left_margin;
  double right = 1.-plot_options_.RightMargin();
  double bottom = 1.-plot_options_.TopMargin();
  double top = 1.;
  if(plot_options_.Title() == TitleType::variable){
    out.push_back(make_shared<TLatex>(0.5*(left+right), 0.5*(bottom+top),
                                      definition_.Title().c_str()));
    out.back()->SetNDC();
    out.back()->SetTextAlign(22);
    out.back()->SetTextFont(plot_options_.Font());
  }else{
    string extra;
    switch(plot_options_.Title()){
    case TitleType::preliminary: extra = "Preliminary"; break;
    case TitleType::simulation: extra = "Simulation"; break;
    case TitleType::supplementary: extra = "Supplementary"; break;
    case TitleType::data: extra = "Data"; break;
    case TitleType::variable:
    default:
      ERROR("Did not understand title type "+to_string(static_cast<int>(plot_options_.Title())));
    }
    out.push_back(make_shared<TLatex>(left, 0.5*(bottom+top),
                                      ("#font[62]{CMS}#scale[0.76]{#font[52]{ "+extra+"}}").c_str()));
    out.back()->SetNDC();
    out.back()->SetTextAlign(12);
    out.back()->SetTextFont(plot_options_.Font());

    ostringstream oss;
    oss << luminosity << " fb^{-1} (13 TeV)" << flush;
    out.push_back(make_shared<TLatex>(right, 0.5*(bottom+top),
                                      oss.str().c_str()));
    out.back()->SetNDC();
    out.back()->SetTextAlign(32);
    out.back()->SetTextFont(plot_options_.Font());
  }
  return out;
}

double HistoStack::FixYAxis(vector<TH1D> &bottom_plots, TPad *top, TPad *bottom){
  double offset = plot_options_.YTitleOffset();
  double margin = plot_options_.LeftMargin();
  if(plot_options_.YAxis() == YAxisType::log){
    margin = 0.15;
    offset = 1.25;
  }else{
    double the_max = GetMaxDraw()*GetLegendRatio();
    int digits = fabs(floor(log10(the_max))-1)+2;

    switch(digits){
    case 0: //Shouldn't happen, but just in case...
    case 1: //Shouldn't happen, but just in case...
    case 2: margin = 0.15; offset = 1.; break;
    default: //For huge axis scale, reverts to scientific notation
    case 3: margin = 0.17; offset = 1.25; break;
    case 4: margin = 0.19; offset = 1.5; break;
    case 5: margin = 0.21; offset = 1.75; break;
    }
  }

  for(auto &h: backgrounds_){
    h.scaled_hist_.SetTitleOffset(offset, "y");
  }
  for(auto &h: signals_){
    h.scaled_hist_.SetTitleOffset(offset, "y");
  }
  for(auto &h: datas_){
    h.scaled_hist_.SetTitleOffset(offset, "y");
  }
  for(auto &h: bottom_plots){
    h.SetTitleOffset(offset, "y");
  }
  top->SetLeftMargin(margin);
  bottom->SetLeftMargin(margin);

  return margin;
}

TLine HistoStack::GetBottomHorizontal() const{
  double left = definition_.Bins().front();
  double right = definition_.Bins().back();
  double y;
  switch(plot_options_.Bottom()){
  case BottomType::ratio: y = 1.; break;
  case BottomType::diff: y = 0.; break;
  case BottomType::off: y = 0.; break;
  default:
    y = 0.;
    DBG("Invalid BottomType: " << to_string(static_cast<int>(plot_options_.Bottom())));
  }

  TLine line(left, y, right, y);
  line.SetNDC(false);
  line.SetLineStyle(3);
  line.SetLineColor(kBlack);
  line.SetLineWidth(2);
  return line;
}

vector<TLine> HistoStack::GetCutLines() const{
  double the_max = GetMaxDraw();
  double the_min = GetMinDraw();

  double bottom;
  switch(plot_options_.YAxis()){
  default:
    DBG("Bad YAxis type " << static_cast<int>(plot_options_.YAxis()));
  case YAxisType::linear: bottom = the_min >= 0. ? 0. : the_min; break;
  case YAxisType::log:    bottom = the_min > 0. ? the_min : plot_options_.LogMinimum(); break;
  }

  vector<TLine> out(definition_.cut_vals_.size());
  for(double cut: definition_.cut_vals_){
    out.emplace_back(cut, bottom, cut, the_max);
    out.back().SetNDC(false);
    out.back().SetLineStyle(2);
    out.back().SetLineColor(kGray);
    out.back().SetLineWidth(3);
  }

  return out;
}

void HistoStack::PrintPlot(double luminosity){
  //Takes already filled histograms and prints to file
  RefreshScaledHistos(luminosity);

  vector<TH1D> bot_plots = GetBottomPlots();
  //I don't know why I can't make this in GetBottomPlots...
  TGraphAsymmErrors bottom_background;
  if(plot_options_.Bottom() != BottomType::off){
    bottom_background = TGraphAsymmErrors(&bot_plots.back());
    bottom_background.SetMinimum(0.1);
    bottom_background.SetMaximum(1.9);
    bot_plots.pop_back();
  }

  TGraphAsymmErrors bkg_error = GetBackgroundError();

  StripTopPlotLabels();
  TLine horizontal = GetBottomHorizontal();
  vector<TLine> cut_vals = GetCutLines();

  unique_ptr<TCanvas> full;
  unique_ptr<TPad> top, bottom;
  GetPads(full, top, bottom);

  double left_margin = plot_options_.LeftMargin();
  if(plot_options_.AutoYAxis()){
    left_margin = FixYAxis(bot_plots, top.get(), bottom.get());
  }

  if(plot_options_.Bottom() != BottomType::off){
    bottom->cd();

    string draw_opt = "ep";
    for(auto &h: bot_plots){
      h.Draw(draw_opt.c_str());
      draw_opt = "ep same";
    }
    bottom_background.Draw("2 same");

    horizontal.Draw("same");

    bottom->RedrawAxis();
    bottom->RedrawAxis("g");
  }

  top->cd();
  if(plot_options_.YAxis() == YAxisType::log) top->SetLogy(true);
  else top->SetLogy(false);

  string draw_opt = "hist";
  DrawAll(backgrounds_, draw_opt);
  if(plot_options_.ShowBackgroundError() && backgrounds_.size()) bkg_error.Draw("2 same");
  DrawAll(signals_, draw_opt);
  ReplaceAll(draw_opt, "hist", "ep");
  DrawAll(datas_, draw_opt);
  for(auto &cut: cut_vals) cut.Draw();

  vector<shared_ptr<TLegend> > legends = GetLegends(left_margin);
  for(auto &legend: legends){
    legend->Draw();
  }

  top->RedrawAxis();
  top->RedrawAxis("g");

  vector<shared_ptr<TLatex> > title_text = GetTitleTexts(luminosity, left_margin);
  for(auto &x: title_text){
    x->Draw();
  }

  string base_name = "plots/"+definition_.Name();
  for(const auto &ext: plot_options_.FileExtensions()){
    string full_name = base_name+"."+ext;
    full->Print(full_name.c_str());
    cout << "Wrote plot to " << full_name << "." << endl;
  }
}

const TH1D & HistoStack::RawHisto(const shared_ptr<Process> &process) const{
  return Histo(process).raw_hist_;
}

TH1D & HistoStack::RawHisto(const shared_ptr<Process> &process){
  return Histo(process).raw_hist_;
}

const TH1D & HistoStack::ScaledHisto(const shared_ptr<Process> &process) const{
  return Histo(process).scaled_hist_;
}

HistoStack & HistoStack::SetPlotOptions(const PlotOpt &plot_opt){
  plot_options_ = plot_opt;
  return *this;
}

const PlotOpt & HistoStack::GetPlotOptions() const{
  return plot_options_;
}

void HistoStack::RefreshScaledHistos(double luminosity){
  StackHistos(luminosity);
  MergeOverflow();
  SetRanges();
  AdjustFillStyles();
}

set<shared_ptr<Process> > HistoStack::GetProcesses() const{
  set<shared_ptr<Process> > processes;
  for(const auto &proc: backgrounds_){
    processes.insert(proc.process_);
  }
  for(const auto &proc: signals_){
    processes.insert(proc.process_);
  }
  for(const auto &proc: datas_){
    processes.insert(proc.process_);
  }
  return processes;
}

std::vector<TH1D> HistoStack::GetBottomPlots() const{
  if(backgrounds_.size() == 0 || plot_options_.Bottom() == BottomType::off){
    return vector<TH1D>();
  }

  TH1D denom = backgrounds_.back().scaled_hist_;
  for(int bin = 0; bin <= denom.GetNbinsX()+1; ++bin){
    denom.SetBinError(bin, 0.);
  }

  vector<TH1D> out(datas_.size()+1);

  for(size_t i = 0; i < datas_.size(); ++i){
    out.at(i) = TH1D(datas_.at(i).scaled_hist_);
  }
  out.back() = TH1D(backgrounds_.back().scaled_hist_);
  out.back().SetFillStyle(3003);
  out.back().SetFillColor(kBlack);
  out.back().SetLineWidth(0);
  out.back().SetMarkerStyle(0);
  out.back().SetMarkerSize(0);

  switch(plot_options_.Bottom()){
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
    ERROR("Bad type for bottom plot: "+to_string(static_cast<int>(plot_options_.Bottom())));
    break;
  }

  double the_min = numeric_limits<double>::infinity();
  double the_max = -numeric_limits<double>::infinity();
  for(auto &h: out){
    h.SetNdivisions(plot_options_.NDivisionsBottom(), "y");
    for(int bin = 1; bin <= h.GetNbinsX(); ++bin){
      double hi = h.GetBinContent(bin)+h.GetBinErrorUp(bin);
      double lo = h.GetBinContent(bin)-fabs(h.GetBinErrorLow(bin));
      if(hi>the_max) the_max = hi;
      if(lo<the_min) the_min = lo;
    }
  }

  if(plot_options_.Bottom() == BottomType::ratio){
    for(auto &h: out){
      h.GetYaxis()->SetTitle("Data/MC");
      h.SetMinimum(0.1);
      h.SetMaximum(1.9);
    }
  }else if(plot_options_.Bottom() == BottomType::diff){
    for(auto &h: out){
      h.GetYaxis()->SetTitle("Data-MC");
      h.SetMinimum(the_min);
      h.SetMaximum(the_max);
    }
  }

  return out;
}

vector<shared_ptr<TLegend> > HistoStack::GetLegends(double left_margin){
  size_t num_procs = backgrounds_.size()+signals_.size()+datas_.size();

  double left = left_margin+plot_options_.LegendPad();
  double right = 1.-plot_options_.RightMargin()-plot_options_.LegendPad();
  double top = 1.-plot_options_.TopMargin()-plot_options_.LegendPad();
  double bottom = top-plot_options_.LegendHeight(num_procs);

  size_t n_entries = datas_.size() + signals_.size() + backgrounds_.size();
  size_t n_columns = min(n_entries, static_cast<size_t>(plot_options_.LegendColumns()));

  double delta_x = (right-left)/n_columns;
  vector<shared_ptr<TLegend> > legends(n_columns);
  for(size_t i = 0; i < n_columns; ++i){
    double x = left+i*delta_x;
    legends.at(i) = make_shared<TLegend>(x, bottom, x+0.12, top);
    legends.at(i)->SetFillStyle(0);
    legends.at(i)->SetBorderSize(0);
    legends.at(i)->SetTextSize(plot_options_.LegendEntryHeight());
    legends.at(i)->SetTextFont(plot_options_.Font());
  }

  size_t entries_added = 0;
  AddEntries(legends, datas_, "lep", n_entries, entries_added);
  AddEntries(legends, signals_, "l", n_entries, entries_added);
  AddEntries(legends, backgrounds_, plot_options_.BackgroundsStacked() ? "f" : "l", n_entries, entries_added);

  return legends;
}

HistoStack::SingleHist::SingleHist(const shared_ptr<Process> &process,
                                   const TH1D &hist):
  process_(process),
  raw_hist_(hist),
  scaled_hist_(){
  raw_hist_.Sumw2();
  scaled_hist_.Sumw2();
}

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

void HistoStack::StackHistos(double luminosity){
  //Scale to luminosity
  for(auto &hist: backgrounds_){
    hist.scaled_hist_ = hist.raw_hist_;
    Scale(hist.scaled_hist_, true);
    hist.scaled_hist_.Scale(luminosity);
  }
  for(auto &hist: signals_){
    hist.scaled_hist_ = hist.raw_hist_;
    Scale(hist.scaled_hist_, true);
    hist.scaled_hist_.Scale(luminosity);
  }
  for(auto &hist: datas_){
    hist.scaled_hist_ = hist.raw_hist_;
    Scale(hist.scaled_hist_, true);
  }

  //Stack histograms
  switch(plot_options_.Stack()){
  case StackType::signal_on_top:
    for(size_t ibkg = 1; ibkg < backgrounds_.size(); ++ibkg){
      backgrounds_.at(ibkg).scaled_hist_ = backgrounds_.at(ibkg).scaled_hist_ + backgrounds_.at(ibkg-1).scaled_hist_;
    }
    if(backgrounds_.size()){
      for(auto &hist: signals_){
        hist.scaled_hist_ = hist.scaled_hist_ + backgrounds_.back().scaled_hist_;
      }
    }
    break;
  case StackType::signal_overlay:
    for(size_t ibkg = 1; ibkg < backgrounds_.size(); ++ibkg){
      backgrounds_.at(ibkg).scaled_hist_ = backgrounds_.at(ibkg).scaled_hist_ + backgrounds_.at(ibkg-1).scaled_hist_;
    }
    break;
  case StackType::lumi_shapes:
    //Don't need to do anything further
    break;
  case StackType::shapes:
    for(auto &hist: backgrounds_){
      Scale(hist.scaled_hist_, false, 100.);
    }
    for(auto &hist: signals_){
      hist.scaled_hist_ = hist.raw_hist_;
      Scale(hist.scaled_hist_, false, 100.);
    }
    for(auto &hist: datas_){
      hist.scaled_hist_ = hist.raw_hist_;
      Scale(hist.scaled_hist_, false, 100.);
    }
    break;
  default:
    DBG("Unknown StackType "+to_string(static_cast<int>(plot_options_.Stack()))+".");
    break;
  }
}

void HistoStack::MergeOverflow(){
  bool underflow = false, overflow = false;
  switch(plot_options_.Overflow()){
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
    ::MergeOverflow(hist.scaled_hist_, underflow, overflow);
  }
  for(auto &hist: signals_){
    ::MergeOverflow(hist.scaled_hist_, underflow, overflow);
  }
  for(auto &hist: datas_){
    ::MergeOverflow(hist.scaled_hist_, underflow, overflow);
  }
}

double HistoStack::GetLegendRatio() const{
  size_t num_plots = backgrounds_.size() + signals_.size() + datas_.size();
  double legend_height = plot_options_.LegendHeight(num_plots);
  double top_plot_height;
  if(plot_options_.Bottom() == BottomType::off){
    top_plot_height = 1.-plot_options_.TopMargin()-plot_options_.BottomMargin();
  }else{
    top_plot_height = 1.-plot_options_.TopMargin()-plot_options_.BottomHeight();
  }
  return top_plot_height/(top_plot_height-legend_height-2.*plot_options_.LegendPad());
}

void HistoStack::SetRanges(){
  double the_min = GetMinDraw();
  double the_max = GetMaxDraw();

  double ratio = GetLegendRatio();

  double top, bottom;
  switch(plot_options_.YAxis()){
  default:
  case YAxisType::linear:
    bottom = the_min >= 0. ? 0. : the_min;
    top = bottom+ratio*(the_max-bottom);
    break;
  case YAxisType::log:
    bottom = the_min > 0. ? the_min : plot_options_.LogMinimum();
    top = exp(log(bottom)+ratio*(log(the_max)-log(bottom)));
    break;
  }

  for(auto &hist: backgrounds_){
    hist.scaled_hist_.SetMinimum(bottom);
    hist.scaled_hist_.SetMaximum(top);
  }
  for(auto &hist: signals_){
    hist.scaled_hist_.SetMinimum(bottom);
    hist.scaled_hist_.SetMaximum(top);
  }
  for(auto &hist: datas_){
    hist.scaled_hist_.SetMinimum(bottom);
    hist.scaled_hist_.SetMaximum(top);
  }
}

void HistoStack::AdjustFillStyles(){
  if(plot_options_.BackgroundsStacked()) return;

  for(auto &bkg: backgrounds_){
    TH1D &h = bkg.scaled_hist_;
    h.SetFillStyle(0);
    h.SetLineColor(h.GetLineColor());
    h.SetLineWidth(5);
  }
}

const vector<HistoStack::SingleHist> & HistoStack::GetHistoList(const shared_ptr<Process> &process) const{
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

vector<HistoStack::SingleHist> & HistoStack::GetHistoList(const shared_ptr<Process> &process){
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

const HistoStack::SingleHist & HistoStack::Histo(const shared_ptr<Process> &process) const{
  const auto &hist_list = GetHistoList(process);
  for(const auto &hist: hist_list){
    if(hist.process_ == process){
      return hist;
    }
  }
  ERROR("Could not find histogram for process "+process->name_+".");
  return hist_list.front();
}

HistoStack::SingleHist & HistoStack::Histo(const shared_ptr<Process> &process){
  auto &hist_list = GetHistoList(process);
  for(auto &hist: hist_list){
    if(hist.process_ == process){
      return hist;
    }
  }
  ERROR("Could not find histogram for process "+process->name_+".");
  return hist_list.front();
}

double HistoStack::GetMaxDraw(double max_bound) const{
  double the_max = -numeric_limits<double>::infinity();
  for(const auto &hist: backgrounds_){
    double this_max = hist.GetMax(max_bound, false);
    if(this_max > the_max && this_max < max_bound){
      the_max = this_max;
    }
  }
  for(const auto &hist: signals_){
    double this_max = hist.GetMax(max_bound, false);
    if(this_max > the_max && this_max < max_bound){
      the_max = this_max;
    }
  }
  for(const auto &hist: datas_){
    double this_max = hist.GetMax(max_bound, true);
    if(this_max > the_max && this_max < max_bound){
      the_max = this_max;
    }
  }
  return the_max;
}

double HistoStack::GetMinDraw(double min_bound) const{
  double the_min = numeric_limits<double>::infinity();
  for(const auto &hist: backgrounds_){
    double this_min = hist.GetMin(min_bound, false);
    if(this_min < the_min && this_min > min_bound){
      the_min = this_min;
    }
  }
  for(const auto &hist: signals_){
    double this_min = hist.GetMin(min_bound, false);
    if(this_min < the_min && this_min > min_bound){
      the_min = this_min;
    }
  }
  for(const auto &hist: datas_){
    double this_min = hist.GetMin(min_bound, true);
    if(this_min < the_min && this_min > min_bound){
      the_min = this_min;
    }
  }
  return the_min;
}

void HistoStack::AddEntries(vector<shared_ptr<TLegend> > &legends,
                            const vector<HistoStack::SingleHist> &hists,
                            const string &style,
                            size_t n_entries,
                            size_t &entries_added) const{
  for(auto h = hists.cbegin(); h != hists.cend(); ++h){
    size_t legend_index = GetLegendIndex(entries_added, n_entries, legends.size());
    TLegend &legend = *legends.at(legend_index);
    string label = h->process_->name_.c_str();
    if(plot_options_.Title() == TitleType::variable){
      double value;
      switch(plot_options_.Stack()){
      default:
        DBG("Bad stack option: " << static_cast<int>(plot_options_.Stack()));
      case StackType::signal_overlay:
      case StackType::signal_on_top:
        value = GetYield(h);
        if(value>=1.){
          label = label + " [N=" + FixedDigits(value, 2) + "]";
        }else{
          label = label + " [N=" + FixedDigits(value, 1) + "]";
        }
        break;
      case StackType::lumi_shapes:
      case StackType::shapes:
        value = GetMean(h);
        label = label + " [#mu=" + FixedDigits(value,3) + "]";
        break;
      }
    }

    legend.AddEntry(&h->scaled_hist_, label.c_str(), style.c_str());
    ++entries_added;
  }
}

double HistoStack::GetYield(vector<HistoStack::SingleHist>::const_iterator h) const{
  TH1D hist = h->scaled_hist_;

  //Subtract underlying histogram
  if(h->process_->type_ == Process::Type::background
     && h != backgrounds_.cbegin()
     && plot_options_.BackgroundsStacked()){
    hist = hist - (--h)->scaled_hist_;
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

double HistoStack::GetMean(vector<HistoStack::SingleHist>::const_iterator h) const{
  TH1D hist = h->scaled_hist_;

  //Subtract underlying histogram
  if(h->process_->type_ == Process::Type::background
     && h != backgrounds_.cbegin()
     && plot_options_.BackgroundsStacked()){
    hist = hist - (--h)->scaled_hist_;
  }

  return hist.GetMean();
}

TGraphAsymmErrors HistoStack::GetBackgroundError() const{
  TGraphAsymmErrors g;
  if(backgrounds_.size() == 0){
    TH1D h("", "", definition_.Nbins(), &definition_.Bins().at(0));
    g = TGraphAsymmErrors(&h);
  }else{
    g = TGraphAsymmErrors(&(backgrounds_.back().scaled_hist_));
  }
  g.SetFillStyle(3003);
  g.SetFillColor(kBlack);
  g.SetLineWidth(0);
  g.SetMarkerSize(0);
  return g;
}
