/*! \class HistoDef

  \brief Contains all information necessary to contruct and fill a histogram,
  independent of drawing style.

  HistoDef provides a means of concisely specifying to other functions and
  classes (for now, mainly HistoStack) the instructions for creating a
  histogram. It contains the binning, the variable to plot, the cut determing
  when to fill the histogram, the weight with which to fill the histogram, the
  x-axis title with units (HistoStack automatically determines the y-axis title,
  but support can be added for manual specification later), and a list of
  significant values of the plotted variable (used by HistoStack to plot
  vertical lines).
 */
#include "histo_def.hpp"

#include <algorithm>

using namespace std;

/*!\brief Constructor using list of bin edges

  \param[in] bins List of Nbins+1 bin edges, including low edge of lowest bin
  and high edge of highest bin

  \param[in] var %Variable with which to fill histogram

  \param[in] x_title X-axis title with format "variable plotted [units]"

  \param[in] cut Cut determining whether histogram is filled

  \param[in] weight Weight with which histogram is filled

  \param[in] cut_vals Values for which to plot a vertical line
*/
HistoDef::HistoDef(const vector<double> &bins,
                   const NamedFunc &var,
                   const string &x_title,
                   const NamedFunc &cut,
                   const NamedFunc &weight,
                   const set<double> &cut_vals):
  tag_(""),
  var_(var),
  cut_(cut),
  weight_(weight),
  x_title_(x_title),
  units_(""),
  cut_vals_(cut_vals),
  bins_(bins){
  sort(bins_.begin(), bins_.end());
  ParseUnits();
  }

/*!\brief Constructor using number of bins and low and high edges of x-axis

  Constructs a histogram with nbins bins linearly spaced from xmin to xmax

  \param[in] nbins Number of bins

  \param[in] xmin Low edge of x-axis

  \param[in] xmax High edge of x-axis

  \param[in] var %Variable with which to fill histogram

  \param[in] x_title X-axis title with format "variable plotted [units]"

  \param[in] cut Cut determining whether histogram is filled

  \param[in] weight Weight with which histogram is filled

  \param[in] cut_vals Values for which to plot a vertical line
*/
HistoDef::HistoDef(size_t nbins,
                   double xmin,
                   double xmax,
                   const NamedFunc &var,
                   const string &x_title,
                   const NamedFunc &cut,
                   const NamedFunc &weight,
                   const set<double> &cut_vals):
  tag_(""),
  var_(var),
  cut_(cut),
  weight_(weight),
  x_title_(x_title),
  units_(""),
  cut_vals_(cut_vals),
  bins_(GetEdges(nbins, xmin, xmax)){
  ParseUnits();
}

/*!\brief Constructor with tag and using list of bin edges

  \param[in] tag Tag for uniquely identifying plot

  \param[in] bins List of Nbins+1 bin edges, including low edge of lowest bin
  and high edge of highest bin

  \param[in] var %Variable with which to fill histogram

  \param[in] x_title X-axis title with format "variable plotted [units]"

  \param[in] cut Cut determining whether histogram is filled

  \param[in] weight Weight with which histogram is filled

  \param[in] cut_vals Values for which to plot a vertical line
*/
HistoDef::HistoDef(const std::string &tag,
		   const vector<double> &bins,
                   const NamedFunc &var,
                   const string &x_title,
                   const NamedFunc &cut,
                   const NamedFunc &weight,
                   const set<double> &cut_vals):
  tag_(tag),
  var_(var),
  cut_(cut),
  weight_(weight),
  x_title_(x_title),
  units_(""),
  cut_vals_(cut_vals),
  bins_(bins){
  sort(bins_.begin(), bins_.end());
  ParseUnits();
  }

/*!\brief Constructor with tag and using number of bins and low and high edges of x-axis

  Constructs a histogram with nbins bins linearly spaced from xmin to xmax

  \param[in] tag Tag for uniquely identifying plot

  \param[in] nbins Number of bins

  \param[in] xmin Low edge of x-axis

  \param[in] xmax High edge of x-axis

  \param[in] var %Variable with which to fill histogram

  \param[in] x_title X-axis title with format "variable plotted [units]"

  \param[in] cut Cut determining whether histogram is filled

  \param[in] weight Weight with which histogram is filled

  \param[in] cut_vals Values for which to plot a vertical line
*/
HistoDef::HistoDef(const std::string &tag,
		   size_t nbins,
                   double xmin,
                   double xmax,
                   const NamedFunc &var,
                   const string &x_title,
                   const NamedFunc &cut,
                   const NamedFunc &weight,
                   const set<double> &cut_vals):
  tag_(tag),
  var_(var),
  cut_(cut),
  weight_(weight),
  x_title_(x_title),
  units_(""),
  cut_vals_(cut_vals),
  bins_(GetEdges(nbins, xmin, xmax)){
  ParseUnits();
}

/*!\brief Get number of bins

  \return Number of bins
*/
size_t HistoDef::Nbins() const{
  return bins_.size()-1;
}

/*!\brief Set binning

  \param[in] bins List of Nbins+1 bin edges, including low edge of lowest bin
  and high edge of highest bin

  \return Reference to *this
*/
HistoDef & HistoDef::Bins(const std::vector<double> &bins){
  bins_ = bins;
  sort(bins_.begin(), bins_.end());
  return *this;
}

/*!\brief Get list of bin edges

  \return List of Nbins+1 bin edges, including low edge of lowest bin and high
  edge of highest bin
*/
const vector<double> & HistoDef::Bins() const{
  return bins_;
}

/*!\brief Get name of plot suitable for use in file name

  \return Output name specifying variable, cut, and weight
*/
string HistoDef::Name() const{
  if(tag_ == ""){
    return var_.PlainName() + "_CUT_" + cut_.PlainName() + "_WGT_" + weight_.PlainName();
  }else{
    return tag_+"_VAR_"+var_.PlainName() + "_CUT_" + cut_.PlainName() + "_WGT_" + weight_.PlainName();
  }
}

/*!\brief Get title of plot containing cut and weight

  Cut is omitted if set to "" or "1." Weight is omitted if set to "weight."

  \return TLatex formatted title with cut and weight
*/
string HistoDef::Title() const{
  bool cut = (cut_.Name() != "" && cut_.Name() != "1");
  bool weight = weight_.Name() != "weight";
  if(cut && weight){
    return cut_.PrettyName()+" (weight="+weight_.PrettyName()+")";
  }else if(cut){
    return cut_.PrettyName();
  }else if(weight){
    return "weight="+weight_.PrettyName();
  }else{
    return "";
  }
}

/*!\brief Create linearly spaced set of bin edges

  \param[in] nbins Number of bins to create

  \param[in] xmin Low edge of lowest bin

  \param[in] xmax High edge of highest bin

  \return List of nbins+1 bin edges linearly spaced from xmin to xmax
*/
vector<double> HistoDef::GetEdges(size_t nbins, double xmin, double xmax){
  vector<double> edges(nbins+1);
  if(nbins != 0){
    double delta = (xmax-xmin)/nbins;
    for(size_t i = 0; i < nbins+1; ++i){
      edges.at(i) = i*delta + xmin;
    }
  }

  //Not necessary, but make sure that first and last edge are correct to available precision
  edges.front() = xmin;
  edges.back() = xmax;

  return edges;
}

/*!\brief Moves units from HistoDef::x_title_ to HistoDef::units_

  Units are found by searching for last occurrence of "[*]"
*/
void HistoDef::ParseUnits(){
  auto p1 = x_title_.rfind('[');
  auto p2 = x_title_.rfind(']');
  if(p1 >= p2 || p2 == string::npos) return;
  units_ = x_title_.substr(p1+1, p2-p1-1);
  x_title_ = x_title_.substr(0, p1);
  while(x_title_.back() == ' '){
    x_title_.pop_back();
  }
}
