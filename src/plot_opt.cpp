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
  title_type_(TitleType::info),
  stack_type_(StackType::signal_overlay),
  overflow_type_(OverflowType::both),
  file_extensions_({"pdf"}),
  title_size_(0.06),
  label_size_(0.045),
  x_title_offset_(1.),
  y_title_offset_(2.2),
  auto_y_axis_(true),
  canvas_width_(600),
  canvas_height_(600),
  left_margin_(0.18),
  right_margin_(0.055),
  bottom_margin_(0.14),
  top_margin_(0.07),
  bottom_height_(0.35),
  legend_columns_(2),
  legend_entry_height_(0.05),
  legend_max_height_(0.3),
  legend_pad_(0.025),
  log_minimum_(0.01),
  n_divisions_(606),
  n_divisions_bottom_(606),
  font_(42),
  show_background_error_(true){
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

PlotOpt & PlotOpt::TitleSize(double title_size){
  title_size_ = title_size;
  return *this;
}

double PlotOpt::TitleSize() const{
  return title_size_;
}

PlotOpt & PlotOpt::LabelSize(double label_size){
  label_size_ = label_size;
  return *this;
}

double PlotOpt::LabelSize() const{
  return label_size_;
}

PlotOpt & PlotOpt::XTitleOffset(double x_title_offset){
  x_title_offset_ = x_title_offset;
  return *this;
}

double PlotOpt::XTitleOffset() const{
  return x_title_offset_;
}

PlotOpt & PlotOpt::YTitleOffset(double y_title_offset){
  y_title_offset_ = y_title_offset;
  return *this;
}

double PlotOpt::YTitleOffset() const{
  return y_title_offset_;
}

PlotOpt & PlotOpt::AutoYAxis(bool auto_y_axis){
  auto_y_axis_ = auto_y_axis;
  return *this;
}

bool PlotOpt::AutoYAxis() const{
  return auto_y_axis_;
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

PlotOpt & PlotOpt::LegendColumns(int columns){
  legend_columns_ = columns;
  return *this;
}

int PlotOpt::LegendColumns() const{
  return legend_columns_;
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

PlotOpt & PlotOpt::LegendPad(double pad){
  legend_pad_ = pad;
  return *this;
}

double PlotOpt::LegendPad() const{
  return legend_pad_;
}

PlotOpt & PlotOpt::NDivisions(int n_divisions){
  n_divisions_ = n_divisions;
  return *this;
}

int PlotOpt::NDivisions() const{
  return n_divisions_;
}

PlotOpt & PlotOpt::NDivisionsBottom(int n_divisions){
  n_divisions_bottom_ = n_divisions;
  return *this;
}

int PlotOpt::NDivisionsBottom() const{
  return n_divisions_bottom_;
}

PlotOpt & PlotOpt::LogMinimum(double log_minimum){
  log_minimum_ = log_minimum;
  return *this;
}

double PlotOpt::LogMinimum() const{
  return log_minimum_;
}

PlotOpt & PlotOpt::Font(int font){
  font_ = font;
  return *this;
}

int PlotOpt::Font() const{
  return font_;
}

PlotOpt & PlotOpt::ShowBackgroundError(bool show_background_error){
  show_background_error_ = show_background_error;
  return *this;
}

bool PlotOpt::ShowBackgroundError() const{
  return show_background_error_;
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

bool PlotOpt::BackgroundsStacked() const{
  switch(stack_type_){
  default:
    DBG("Invalid stack type " << static_cast<int>(stack_type_));
  case StackType::signal_overlay:
  case StackType::signal_on_top:
    return true;
  case StackType::lumi_shapes:
  case StackType::shapes:
    return false;
  }
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
  }else if(property == "FileExtensions"){
    auto exts = FileExtensions();
    exts.insert(value);
    FileExtensions(exts);
  }else if(property == "TitleSize"){
    TitleSize(stod(value));
  }else if(property == "LabelSize"){
    LabelSize(stod(value));
  }else if(property == "xTitleOffset" || property == "XTitleOffset"){
    XTitleOffset(stod(value));
  }else if(property == "yTitleOffset" || property == "YTitleOffset"){
    YTitleOffset(stod(value));
  }else if(property == "AutoYAxis"){
    AutoYAxis(static_cast<bool>(stoi(value)));
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
  }else if(property == "LegendColumns"){
    LegendColumns(stoi(value));
  }else if(property == "LegendEntrySize" || property == "LegendSize"){
    LegendEntryHeight(stod(value));
  }else if(property == "LegendMaxSize"){
    LegendMaxHeight(stod(value));
  }else if(property == "LegendPad"){
    LegendPad(stod(value));
  }else if(property == "BottomPlotHeight"){
    BottomHeight(stod(value));
  }else if(property == "LogMinimum"){
    LogMinimum(stod(value));
  }else if(property == "nDivisions" || property == "NDivisions"){
    NDivisions(stoi(value));
  }else if(property == "nDivisionsBottom" || property == "NDivisionsBottom"){
    NDivisionsBottom(stoi(value));
  }else if(property == "Font"){
    Font(stoi(value));
  }else if(property == "ShowBackgroundError"){
    ShowBackgroundError(stoi(value));
  }else{
    DBG("Did not understand property name "<<property);
  }
}
