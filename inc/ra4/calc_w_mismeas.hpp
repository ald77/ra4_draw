#ifndef H_CALC_W_MISMEAS
#define H_CALC_W_MISMEAS

#include <vector>
#include <string>

#include "core/gamma_params.hpp"
#include "core/named_func.hpp"

void PrintYield(const std::string &name, const std::vector<GammaParams> &yields);

bool UseData();

NamedFunc MismeasurementCut();

NamedFunc MismeasurementWeight();

NamedFunc WeightWithCut();

void BayesianOmega(double &omega_lo, double &omega_mid, double &omega_hi,
                   const GammaParams &a_data, const GammaParams &b_data,
                   const GammaParams &c_data, const GammaParams &d_data,
                   const GammaParams &a_mc_good, const GammaParams &b_mc_good,
                   const GammaParams &c_mc_good, const GammaParams &d_mc_good,
                   const GammaParams &a_mc_bad, const GammaParams &b_mc_bad,
                   const GammaParams &c_mc_bad, const GammaParams &d_mc_bad);

double GetOmega(double mua_data, double mub_data, double muc_data, double mud_data,
                double mua_mc_good, double mub_mc_good, double muc_mc_good, double mud_mc_good,
                double mua_mc_bad, double mub_mc_bad, double muc_mc_bad, double mud_mc_bad);

void GetInterval(double &x_lo, double &x_hi, const std::vector<double> &x_vals);

void GetOptions(int argc, char *argv[]);

#endif
