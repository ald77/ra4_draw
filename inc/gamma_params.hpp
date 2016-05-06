#ifndef H_GAMMA_PARAM
#define H_GAMMA_PARAM

#include <ostream>

class GammaParams{
public:
  GammaParams();
  GammaParams(double n_effective, double weight);
  GammaParams(const GammaParams &) = default;
  GammaParams & operator=(const GammaParams &) = default;
  GammaParams(GammaParams &&) = default;
  GammaParams & operator=(GammaParams &&) = default;
  ~GammaParams() = default;

  void SetNEffectiveAndWeight(double n_effective, double weight);
  void SetYieldAndUncertainty(double yield, double uncertainty);

  double Yield() const;
  void Yield(double yield);

  double Uncertainty() const;
  void Uncertainty(double uncertainty);

  double NEffective() const;
  void NEffective(double n_effective);

  double Weight() const;
  void Weight(double weight);

  double CorrectedUncertainty() const;

  GammaParams & operator+=(const GammaParams &gp);
  GammaParams & operator*=(double scale);

private:
  double n_effective_, weight_;
};

GammaParams operator+(GammaParams gp1, GammaParams gp2);
GammaParams operator*(double scale, GammaParams gp);
GammaParams operator*(GammaParams gp, double scale);
std::ostream & operator<<(std::ostream &stream, const GammaParams &gp);

#endif
