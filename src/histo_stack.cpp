#include "histo_stack.hpp"

#include <algorithm>
#include <sstream>

#include "TCanvas.h"

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

  void ResizeLabels(TH1D &h, double scale){
    h.GetXaxis()->SetTitle("");
    h.SetLabelSize(0., "x");
    h.SetTitleSize(0., "x");
    return;
    h.SetMarkerSize(h.GetMarkerSize()*scale);
    h.SetLabelSize(h.GetLabelSize("x")*scale, "x");
    h.SetLabelSize(h.GetLabelSize("y")*scale, "y");
    h.SetLabelSize(h.GetLabelSize("z")*scale, "z");
    h.SetTitleSize(h.GetTitleSize("x")*scale, "x");
    h.SetTitleSize(h.GetTitleSize("y")*scale, "y");
    h.SetTitleSize(h.GetTitleSize("z")*scale, "z");
    h.SetLabelOffset(h.GetLabelOffset("x")/scale, "x");
    h.SetLabelOffset(h.GetLabelOffset("y")/scale, "y");
    h.SetLabelOffset(h.GetLabelOffset("z")/scale, "x");
    h.SetTitleOffset(h.GetTitleOffset("x")/scale, "x");
    h.SetTitleOffset(h.GetTitleOffset("y")/scale, "y");
    h.SetTitleOffset(h.GetTitleOffset("z")/scale, "x");
    h.SetTickLength(h.GetTickLength("x")*scale, "x");
    //h.SetTickLength(h.GetTickLength("y")/scale, "y");
    //h.SetTickLength(h.GetTickLength("z")*scale, "z");
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

  double bin_width = (definition.bins_.back()-definition.bins_.front())/(definition.bins_.size()-1);
  
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
  
  TH1D empty("", (";"+x_title+";"+oss.str()).c_str(), definition.GetNbins(), &definition.bins_.at(0));
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

void HistoStack::ResizeTopPlotLabels(){
  if(plot_options_.Bottom() == BottomType::off) return;
  double scale = 1./(1.-plot_options_.BottomHeight());
  for(auto &hist: backgrounds_){
    ResizeLabels(hist.scaled_hist_, scale);
  }
  for(auto &hist: signals_){
    ResizeLabels(hist.scaled_hist_, scale);
  }
  for(auto &hist: datas_){
    ResizeLabels(hist.scaled_hist_, scale);
  }
}

void HistoStack::PrintPlot(){
  //Takes already filled histograms and prints to file
  RefreshScaledHistos();  
  vector<TH1D> bot_plots = GetBottomPlots();
  ResizeTopPlotLabels();
  
  TCanvas full("", "",
               plot_options_.CanvasWidth(),
               plot_options_.CanvasHeight());
  full.cd();
  unique_ptr<TPad> top, bottom;
  full.SetMargin(plot_options_.LeftMargin(),
                 plot_options_.RightMargin(),
                 plot_options_.BottomMargin(),
                 plot_options_.TopMargin());
  full.SetMargin(0., 0., 0., 0.);
  full.SetTicks(1,1);
  full.SetFillStyle(4000);
  if(plot_options_.Bottom() == BottomType::off){
    top.reset(new TPad("", "", 0., 0., 1., 1.));
    bottom.reset();
    top->SetTicks(1,1);
    top->SetFillStyle(4000);
    top->SetMargin(plot_options_.LeftMargin(),
                   plot_options_.RightMargin(),
                   plot_options_.BottomMargin(),
                   plot_options_.TopMargin());
    top->Draw();
  }else{
    top.reset(new TPad("", "", 0., 0., 1., 1.));
    bottom.reset(new TPad("", "", 0., 0., 1., 1.));
    top->SetTicks(1,1);
    bottom->SetTicks(1,1);
    top->SetFillStyle(4000);
    bottom->SetFillStyle(4000);
    top->SetMargin(plot_options_.LeftMargin(),
                   plot_options_.RightMargin(),
                   plot_options_.BottomHeight(),
                   plot_options_.TopMargin());///(1.-plot_options_.BottomHeight()));
    bottom->SetMargin(plot_options_.LeftMargin(),
                      plot_options_.RightMargin(),
                      plot_options_.BottomMargin(),
                      1.-plot_options_.BottomHeight());
    bottom->Draw();
    top->Draw();
  }
  
  string draw_opt = "hist";
  if(bottom){
    bottom->cd();
    
    for(auto &h: bot_plots){
      h.Draw(draw_opt.c_str());
      draw_opt = "hist same";
      break;
    }
    
    bottom->RedrawAxis();
    bottom->RedrawAxis("g");
    draw_opt = "hist";
  }

  top->cd();
  if(plot_options_.YAxis() == YAxisType::log){
    top->SetLogy(true);
  }else{
    top->SetLogy(false);
  }

  DrawAll(backgrounds_, draw_opt);
  DrawAll(signals_, draw_opt);
  ReplaceAll(draw_opt, "hist", "ep");
  DrawAll(datas_, draw_opt);
  auto legend = GetLegend();
  legend->Draw();
  top->RedrawAxis();
  top->RedrawAxis("g");

  string base_name = "plots/"+definition_.GetName();
  for(const auto &ext: plot_options_.FileExtensions()){
    string full_name = base_name+"."+ext;
    full.Print(full_name.c_str());
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

void HistoStack::RefreshScaledHistos(){
  StackHistos();
  MergeOverflow();
  SetRanges();
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
  if(backgrounds_.size() == 0) return vector<TH1D>();

  TH1D denom = backgrounds_.front().scaled_hist_;
  for(int bin = 0; bin <= denom.GetNbinsX()+1; ++bin){
    denom.SetBinError(bin, 0.);
  }

  vector<TH1D> out(datas_.size()+1);
  
  for(size_t i = 0; i < datas_.size(); ++i){
    out.at(i) = TH1D(datas_.at(i).scaled_hist_);
    out.at(i).Divide(&denom);
  }
  out.back() = TH1D(backgrounds_.back().scaled_hist_);
  out.back().Divide(&denom);

  double the_min = numeric_limits<double>::infinity();
  double the_max = -numeric_limits<double>::infinity();
  for(auto &h: out){
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

unique_ptr<TLegend> HistoStack::GetLegend(){
  size_t num_procs = backgrounds_.size()+signals_.size()+datas_.size();
  
  double left = plot_options_.LeftMargin()+plot_options_.LegendPad();
  double right = 1.-plot_options_.RightMargin()-plot_options_.LegendPad();
  double top = 1.-plot_options_.TopMargin()-plot_options_.LegendPad();
  double bottom = top-plot_options_.LegendHeight(num_procs);
  //top = plot_options_.GlobalToTopYNDC(top);
  //bottom = plot_options_.GlobalToTopYNDC(bottom);

  auto legend = unique_ptr<TLegend>(new TLegend(left, bottom, right, top));
  legend->SetNColumns(2);
  legend->SetFillStyle(0);
  legend->SetFillColor(0);
  legend->SetBorderSize(0);
  
  for(const auto &hist: datas_){
    legend->AddEntry(&hist.scaled_hist_, hist.process_->name_.c_str(), "e0p");
  }
  for(const auto &hist: backgrounds_){
    legend->AddEntry(&hist.scaled_hist_, hist.process_->name_.c_str(), "f");
  }
  for(const auto &hist: signals_){
    legend->AddEntry(&hist.scaled_hist_, hist.process_->name_.c_str(), "l");
  }

  return legend;
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

void HistoStack::StackHistos(){
  for(auto &hist: backgrounds_){
    hist.scaled_hist_ = hist.raw_hist_;
    Scale(hist.scaled_hist_, true);
  }
  for(auto &hist: signals_){
    hist.scaled_hist_ = hist.raw_hist_;
    Scale(hist.scaled_hist_, true);
  }
  for(auto &hist: datas_){
    hist.scaled_hist_ = hist.raw_hist_;
    Scale(hist.scaled_hist_, true);
  }
  
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

void HistoStack::SetRanges(){
  double the_min = GetMinDraw();
  double the_max = GetMaxDraw();

  size_t num_plots = backgrounds_.size() + signals_.size() + datas_.size();
  double legend_height = plot_options_.LegendHeight(num_plots);
  double top_plot_height;
  if(plot_options_.Bottom() == BottomType::off){
    top_plot_height = 1.-plot_options_.TopMargin()-plot_options_.BottomMargin();
  }else{
    top_plot_height = 1.-plot_options_.TopMargin()-plot_options_.BottomHeight();
  }
  double ratio = top_plot_height/(top_plot_height-legend_height-2.*plot_options_.LegendPad());

  double bottom, top;
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
