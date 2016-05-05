#include "plot_opt.hpp"

#include <cmath>

#include <algorithm>
#include <fstream>

#include "utilities.hpp"

using namespace std;
using namespace PlotOptTypes;

PlotOpt::PlotOpt():
  bottom_type_(BottomType::off),
  y_axis_type_(YAxisType::linear),
  title_type_(TitleType::variable),
  stack_type_(StackType::signal_overlay),
  overflow_type_(OverflowType::both),
  file_extensions_({"pdf"}),
  canvas_width_(800),
  canvas_height_(800),
  left_margin_(0.15),
  right_margin_(0.05),
  bottom_margin_(0.1),
  top_margin_(0.1),
  bottom_height_(1./3.),
  legend_entry_height_(0.05),
  legend_max_height_(0.3){
}

PlotOpt::PlotOpt(const string &file_name,
                 const string &config_name):
  PlotOpt(){
  LoadOptions(file_name, config_name);
}

PlotOpt PlotOpt::operator()() const{
  return *this;
}

PlotOpt & PlotOpt::LoadOptions(const string &file_name,
                               const string &config_name){
  ifstream file(file_name);
  string line;
  string current_config = "";
  int line_num = 0;
  while(getline(file, line)){
    ++line_num;
    ReplaceAll(line, " ", "");
    ReplaceAll(line, "\t", "");
    auto start  = line.find('[');
    auto end = line.find(']');
    if(start==string::npos && end!=string::npos){
      ERROR("Could not find closing brace in line "+to_string(line_num));
    }
    if(start!=string::npos && end==string::npos){
      ERROR("Could not find closing brace in line "+to_string(line_num));
    }
    if(start<end && start != string::npos && end != string::npos){
      current_config = line.substr(start+1, end-start-1);
    }else if(current_config == config_name
             && line.size()
             && line.at(0)!='#'){
      auto pos = line.find("=");
      if(pos == string::npos) continue;
      string prop_name = line.substr(0,pos);
      string value = line.substr(pos+1);
      SetProperty(prop_name, value);
    }
  }
  return *this;
}

PlotOpt & PlotOpt::Bottom(BottomType bottom_type){
  bottom_type_ = bottom_type;
  return *this;
}

BottomType PlotOpt::Bottom() const{
  return bottom_type_;
}

PlotOpt & PlotOpt::YAxis(YAxisType y_axis_type){
  y_axis_type_ = y_axis_type;
  return *this;
}

YAxisType PlotOpt::YAxis() const{
  return y_axis_type_;  
}

PlotOpt & PlotOpt::Title(TitleType title_type){
  title_type_ = title_type;
  return *this;
}

TitleType PlotOpt::Title() const{
  return title_type_;
}

PlotOpt & PlotOpt::Stack(StackType stack_type){
  stack_type_ = stack_type;
  return *this;
}

StackType PlotOpt::Stack() const{
  return stack_type_;
}

PlotOpt & PlotOpt::Overflow(OverflowType overflow_type){
  overflow_type_ = overflow_type;
  return *this;
}

OverflowType PlotOpt::Overflow() const{
  return overflow_type_;
}

PlotOpt & PlotOpt::FileExtensions(const set<string> &file_extensions){
  file_extensions_ = file_extensions;
  return *this;
}

const set<string> & PlotOpt::FileExtensions() const{
  return file_extensions_;
}

PlotOpt & PlotOpt::CanvasSize(int width, int height){
  canvas_width_ = width;
  canvas_height_ = height;
  return *this;
}

PlotOpt & PlotOpt::CanvasWidth(int width){
  canvas_width_ = width;
  return *this;
}

int PlotOpt::CanvasWidth() const{
  return canvas_width_;
}

PlotOpt & PlotOpt::CanvasHeight(int height){
  canvas_height_ = height;
  return *this;
}

int PlotOpt::CanvasHeight() const{
  return canvas_height_;
}

PlotOpt & PlotOpt::Margin(double left, double right, double bottom, double top){
  left_margin_ = left;
  right_margin_ = right;
  bottom_margin_ = bottom;
  top_margin_ = top;
  return *this;
}

PlotOpt & PlotOpt::LeftMargin(double left){
  left_margin_ = left;
  return *this;
}

double PlotOpt::LeftMargin() const{
  return left_margin_;
}

PlotOpt & PlotOpt::RightMargin(double right){
  right_margin_ = right;
  return *this;
}

double PlotOpt::RightMargin() const{
  return right_margin_;
}

PlotOpt & PlotOpt::BottomMargin(double bottom){
  bottom_margin_ = bottom;
  return *this;
}

double PlotOpt::BottomMargin() const{
  return bottom_margin_;
}

PlotOpt & PlotOpt::TopMargin(double top){
  top_margin_ = top;
  return *this;
}

double PlotOpt::TopMargin() const{
  return top_margin_;
}

PlotOpt & PlotOpt::BottomHeight(double bottom_height){
  bottom_height_ = bottom_height;
  return *this;
}

double PlotOpt::BottomHeight() const{
  return bottom_height_;
}

PlotOpt & PlotOpt::LegendEntryHeight(double height){
  legend_entry_height_ = height;
  return *this;
}

double PlotOpt::LegendEntryHeight() const{
  return legend_entry_height_;
}

PlotOpt & PlotOpt::LegendMaxHeight(double height){
  legend_max_height_ = height;
  return *this;
}

double PlotOpt::LegendMaxHeight() const{
  return legend_max_height_;
}

double PlotOpt::TopToGlobalYNDC(double top_y) const{
  if(bottom_type_ == BottomType::off) return top_y;
  else return 1.-(1.-top_y)*(1. - bottom_margin_ - bottom_height_);
}

double PlotOpt::GlobalToTopYNDC(double global_y) const{
  if(bottom_type_ == BottomType::off) return global_y;
  else return 1.-(1.-global_y)/(1. - bottom_margin_ - bottom_height_);
}

double PlotOpt::BottomToGlobalYNDC(double bottom_y) const{
  if(bottom_type_ == BottomType::off) return bottom_y;
  else return bottom_y*(bottom_margin_+bottom_height_);
}

double PlotOpt::GlobalToBottomYNDC(double global_y) const{
  if(bottom_type_ == BottomType::off) return global_y;
  else return global_y*(bottom_margin_+bottom_height_);
}

double PlotOpt::LegendHeight(size_t num_entries) const{
  return min(legend_max_height_, legend_entry_height_*ceil(num_entries*0.5));
}

void PlotOpt::SetProperty(const string &property,
                          const string &value){
  if(property == "BottomType"){
    Bottom(static_cast<BottomType>(stoi(value)));
  }else if(property == "YAxisType"){
    YAxis(static_cast<YAxisType>(stoi(value)));
  }else if(property == "TitleType"){
    Title(static_cast<TitleType>(stoi(value)));
  }else if(property == "StackType"){
    Stack(static_cast<StackType>(stoi(value)));
  }else if(property == "OverflowType"){
    Overflow(static_cast<OverflowType>(stoi(value)));
  }else if(property == "CanvasWidth" || property == "CanvasW"){
    CanvasWidth(stoi(value));
  }else if(property == "CanvasHeight" || property == "CanvasH"){
    CanvasHeight(stoi(value));
  }else if(property == "PadLeftMargin"){
    LeftMargin(stod(value));
  }else if(property == "PadRightMargin"){
    RightMargin(stod(value));
  }else if(property == "PadBottomMargin"){
    BottomMargin(stod(value));
  }else if(property == "PadTopMargin"){
    TopMargin(stod(value));
  }else if(property == "LegendEntrySize" || property == "LegendSize"){
    LegendEntryHeight(stod(value));
  }else if(property == "LegendMaxSize"){
    LegendMaxHeight(stod(value));
  }else if(property == "BottomPlotHeight"){
    BottomHeight(stod(value));
  }else{
    DBG("Did not understand property name "<<property);
  }
}
