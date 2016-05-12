#ifndef H_HISTO_DEF
#define H_HISTO_DEF

#include <vector>
#include <set>
#include <string>

#include "named_func.hpp"

class HistoDef{
public:
  HistoDef(const std::vector<double> &bins,
           const NamedFunc &var,
           const std::string &x_title = "",
           const NamedFunc &cut = 1.,
           const NamedFunc &weight = "weight",
           const std::set<double> &cut_vals = {});
  HistoDef(size_t nbins,
           double xmin,
           double xmax,
           const NamedFunc &var,
           const std::string &x_title = "",
           const NamedFunc &cut = 1.,
           const NamedFunc &weight = "weight",
           const std::set<double> &cut_vals = {});
  HistoDef(const HistoDef &) = default;
  HistoDef& operator=(const HistoDef &) = default;
  HistoDef(HistoDef &&) = default;
  HistoDef& operator=(HistoDef &&) = default;
  ~HistoDef() = default;

  size_t Nbins() const;
  HistoDef & Bins(const std::vector<double> &bins);
  const std::vector<double> & Bins() const;

  std::string Name() const;
  std::string Title() const;

  NamedFunc var_;//!< %Variable to be plotted
  NamedFunc cut_;//!< Cut determining whether histogram is filled
  NamedFunc weight_;//!< Weight with which to fill histogram
  std::string x_title_;//!< X-axis title _without_ units
  std::string units_;//!< Units of HistoDef::var_
  std::set<double> cut_vals_;//!< Values of HistoDef::var_ for which to plot a vertical line

private:
  std::vector<double> bins_;
  static std::vector<double> GetEdges(size_t nbins, double xmin, double xmax);
  void ParseUnits();
};

#endif
