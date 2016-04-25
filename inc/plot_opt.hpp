#ifndef H_PLOT_OPT
#define H_PLOT_OPT

#include <cstddef>

#include <string>

namespace PlotOptTypes{
  enum class BottomType{off, ratio, diff, signif};
  enum class YAxisType{linear, log};
  enum class TitleType{variable, preliminary, simulation, supplementary, data};
  enum class StackType{signal_overlay, signal_on_top, lumi_shapes, shapes};
  enum class OverflowType{none, underflow, overflow, both};
}

class PlotOpt{
public:
  PlotOpt();
  PlotOpt(const std::string &file_name,
          const std::string &config_name);
  PlotOpt(const PlotOpt &) = default;
  PlotOpt& operator=(const PlotOpt &) = default;
  PlotOpt(PlotOpt &&) = default;
  PlotOpt& operator=(PlotOpt &&) = default;
  ~PlotOpt() = default;

  PlotOpt operator()() const;

  PlotOpt & LoadOptions(const std::string &file_name,
                        const std::string &config_name);

  PlotOpt & Bottom(PlotOptTypes::BottomType bottom_type);
  PlotOptTypes::BottomType Bottom() const;

  PlotOpt & YAxis(PlotOptTypes::YAxisType y_axis_type);
  PlotOptTypes::YAxisType YAxis() const;

  PlotOpt & Title(PlotOptTypes::TitleType title_type);
  PlotOptTypes::TitleType Title() const;

  PlotOpt & Stack(PlotOptTypes::StackType stack_type);
  PlotOptTypes::StackType Stack() const;

  PlotOpt & Overflow(PlotOptTypes::OverflowType overflow_type);
  PlotOptTypes::OverflowType Overflow() const;

  PlotOpt & CanvasSize(int width, int height);
  PlotOpt & CanvasWidth(int width);
  int CanvasWidth() const;
  PlotOpt & CanvasHeight(int height);
  int CanvasHeight() const;

  PlotOpt & Margin(double left, double right, double bottom, double top);
  PlotOpt & LeftMargin(double left);
  double LeftMargin() const;
  PlotOpt & RightMargin(double right);
  double RightMargin() const;
  PlotOpt & BottomMargin(double bottom);
  double BottomMargin() const;
  PlotOpt & TopMargin(double top);
  double TopMargin() const;

  PlotOpt & BottomHeight(double bottom_height);
  double BottomHeight() const;

  PlotOpt & LegendEntryHeight(double height);
  double LegendEntryHeight() const;
  PlotOpt & LegendMaxHeight(double height);
  double LegendMaxHeight() const;

  double TopToGlobalYNDC(double top_y) const;
  double GlobalToTopYNDC(double global_y) const;
  double BottomToGlobalYNDC(double bottom_y) const;
  double GlobalToBottomYNDC(double global_y) const;

  double LegendHeight(size_t num_entries) const;

private:
  PlotOptTypes::BottomType bottom_type_;
  PlotOptTypes::YAxisType y_axis_type_;
  PlotOptTypes::TitleType title_type_;
  PlotOptTypes::StackType stack_type_;
  PlotOptTypes::OverflowType overflow_type_;
  int canvas_width_, canvas_height_;
  double left_margin_, right_margin_, bottom_margin_, top_margin_;
  double bottom_height_;
  double legend_entry_height_, legend_max_height_;

  void SetProperty(const std::string &property_name,
                   const std::string &value_string);
};

#endif
