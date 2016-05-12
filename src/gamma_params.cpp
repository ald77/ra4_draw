/*! \class GammaParams

  \brief Represents a yield and uncertainty obtained by counting weighted events

  GammaParams uses two different representations of rate: one using the
  rate/yield and uncertainty and another using the number of unweighted events
  and the weight. In the case that all events have the same weight, the
  conversion is:

  - Yield or rate = (unweighted events) * (weight)

  - Uncertainty = sqrt(unweighted events) * weight

  - Unweighted events = (yield / uncertainty)^2

  - Weight = uncertainty^2/yield

  In the case of events having different weights, GammaParams uses an
  approximation to get an effective number of unweighted events and an effective
  weight for the combination in the same way as ROOT's Sumw2;

  The GammaParams::Uncertainty() method returns the "sqrt(N)"-like uncertainty
  on a Poisson mean correspdoning to a prior of Gamma(alpha=0, beta=0). This is
  the uncertainty used by the uncertainty getter and setter methods. The
  GammaParams::CorrectedUncertainty() returns the "sqrt(N+1)"-like uncertainty
  corresponding to a prior of Gamma(alpha=1, beta=0).
*/

#include "gamma_params.hpp"

#include <cmath>

/*!\brief Standard constructor initializing to unweighted count = 0, weight = 0
 */
GammaParams::GammaParams():
  n_effective_(0.),
  weight_(0.){
}

/*!\brief Constructs GammaParams with given number of unweighted events and
  weight

  \param[in] n_effective Effective number of unweighted events

  \param[in] weight Effective weight of sample
*/
GammaParams::GammaParams(double n_effective, double weight):
  n_effective_(n_effective),
  weight_(weight){
  }

/*!\brief Sets GammaParams by yield and uncertainty

  Resets the number of unweighted events and effective weight to corresponding
  values.

  \param[in] yield The _weighted_ number of events

  \param[in] uncertainty The "sqrt(N)"-like uncertainty on the yield/rate.

  \see GammaParams::Uncertainty()
*/
void GammaParams::SetYieldAndUncertainty(double yield, double uncertainty){
  if(yield > 0.){
    n_effective_ = (yield*yield)/(uncertainty*uncertainty);
    weight_ = uncertainty*uncertainty/yield;
  }else{
    n_effective_ = 0.;
    weight_ = hypot(yield, uncertainty);
  }
}

/*!\brief Sets GammaParams by unweighted event count and weight

  Resets the yield/rate and uncertainty to corresponding values.

  \param[in] n_effective The effective number of _unweighted_events. It need not
  be an integer.

  \param[in] weight The effective weight of the sample
*/
void GammaParams::SetNEffectiveAndWeight(double n_effective, double weight){
  n_effective_ = n_effective;
  weight_ = weight;
}

/*!\brief Get the yield/rate

  \return The yield or rate
*/
double GammaParams::Yield() const{
  return n_effective_*weight_;
}

/*!\brief Set the yield/rate

  \param[in] yield The _weighted_ number of events
*/
void GammaParams::Yield(double yield){
  SetYieldAndUncertainty(yield, Uncertainty());
}

/*!\brief Get the "sqrt(N)"-like uncertainty on the yield/rate

  \return The "sqrt(N)"-like uncertainty corresponding to a prior of
  Gamma(alpha=0, beta=0) on the Poisson rate. This is the uncertainty used by
  the GammaParams setter methods.
*/
double GammaParams::Uncertainty() const{
  return sqrt(n_effective_)*weight_;
}

/*!\brief Set the "sqrt(N)"-like uncertainty on the yield/rate

  \param[in] uncertainty The "sqrt(N)"-like uncertainty corresponding to a prior
  of Gamma(alpha=0, beta=0) on the Poisson rate. This is the uncertainty used by
  GammaParams::Uncertainty().
*/
void GammaParams::Uncertainty(double uncertainty){
  SetYieldAndUncertainty(Yield(), uncertainty);
}

/*!\brief Get the effective number of unweighted events

  \return Effective number of unweighted events
*/
double GammaParams::NEffective() const{
  return n_effective_;
}

/*!\brief Set the effective number of unweighted events

  \param[in] n_effective Effective number of unweighted events
*/
void GammaParams::NEffective(double n_effective){
  SetNEffectiveAndWeight(n_effective, Weight());
}

/*!\brief Get the effective weight for the sample

  \return The effective weight for the sample
*/
double GammaParams::Weight() const{
  return weight_;
}

/*!\brief Set the effective weight for the sample

  \param[in] weight The effective weight for the sample
*/
void GammaParams::Weight(double weight){
  SetNEffectiveAndWeight(NEffective(), weight);
}

/*!\brief Get the "sqrt(N+1)"-like uncertainty on the yield/rate

  \return The "sqrt(N+1)"-like uncertainty corresponding to a prior of
  Gamma(alpha=1, beta=0) on the Poisson rate. This is _not_ the uncertainty used
  by the GammaParams setter methods.
*/
double GammaParams::CorrectedUncertainty() const{
  return sqrt(n_effective_+1.)*weight_;
}

/*!\brief Adds another GammaParams to *this

  Uses a Sumw2-like approximation to obtain effective number of unweighted
  events and weight. The resulting yield is the sum of the original two
  yields. The resulting uncertainty is the square-root-of-the-sum-of-squares of
  the original uncertainties.

  \param[in] gp GammaParams to add to *this

  \return Reference to *this
*/
GammaParams & GammaParams::operator+=(const GammaParams &gp){
  if(NEffective() == 0. && gp.NEffective() == 0.){
    SetNEffectiveAndWeight(0., hypot(Weight(), gp.Weight()));
  }else{
    SetYieldAndUncertainty(Yield()+gp.Yield(), hypot(Uncertainty(),gp.Uncertainty()));
  }
  return *this;
}

/*!\brief Scale/multiply *this by a constant

  The yield and uncertainty scale linearly.

  \param[in] scale Amount by which to scale *this

  \return Reference to *this
*/
GammaParams & GammaParams::operator*=(double scale){
  SetNEffectiveAndWeight(NEffective(), scale*Weight());
  return *this;
}

/*!\brief Add two GammaParams

  \param[in] gp1 Left-hand addend

  \param[in] gp2 Right-hand addend

  \return Regerence to *this

  \see GammaParams::operator+=()
*/
GammaParams operator+(GammaParams gp1, GammaParams gp2){
  return (gp1 += gp2);
}

/*!\brief Scale/multiply a GammaParams by a constant on the left

  \param[in] scale Amount by which to scale gp

  \param[in] GammaParams to be scaled

  \return Scaled copy of gp

  \see GammaParams::operator*=()
*/
GammaParams operator*(double scale, GammaParams gp){
  return (gp*=scale);
}

/*!\brief Scale/multiply a GammaParams by a constant on the right

  \param[in] GammaParams to be scaled

  \param[in] scale Amount by which to scale gp

  \return Scaled copy of gp

  \see GammaParams::operator*=()
*/
GammaParams operator*(GammaParams gp, double scale){
  return (gp*=scale);
}

/*!\brief Print GammaParams to output stream

  \param[in,out] stream Output stream to print to

  \param[in] fp GammaParams to print

  \return Reference to stream
*/
std::ostream & operator<<(std::ostream &stream, const GammaParams &gp){
  stream << gp.Yield() << "+-" << gp.CorrectedUncertainty()
         << " (N=" << gp.NEffective() << ", w=" << gp.Weight() << ")";
  return stream;
}
