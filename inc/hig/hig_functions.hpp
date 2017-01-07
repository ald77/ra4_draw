#ifndef H_HIG_FUNCTIONS
#define H_HIG_FUNCTIONS

#include <cstddef>

#include <string>

#include "core/named_func.hpp"
 
namespace Higfuncs{
	// count number of jets associated with true B-hadrons 
	// among the 4 jets used to build average mjj
	extern const NamedFunc ntrub;

	// nominal and control region b-tag categorization
	extern const NamedFunc hig_nb;
	extern const NamedFunc hig_nb_extended;

	// alternative b-tag category options
	extern const NamedFunc hig_nb_mmmm;
	extern const NamedFunc hig_nb_ttll;
	extern const NamedFunc hig_nb_tmml;

	// weights derived in data/MC comparisons of MET and nb in CRs
	NamedFunc::ScalarType wgt_nb_met(const Baby &b, bool ttonly);
	// weight that allows subtracting the reweighted ttbar from the data histogram
	NamedFunc::ScalarType wgt_subtr_ttx(const Baby &b, std::string json);
	// namedfunc to apply the weights in wgt_nb_met to all bkg. processes
	extern const NamedFunc wgt_comp;

	// calculate effect of systematics calculated for each background 
	// in the data control regions on the total bkg. kappa
	extern const NamedFunc wgt_syst_ttx;
	extern const NamedFunc wgt_syst_vjets;
	extern const NamedFunc wgt_syst_qcd;


	// analysis trigger and its efficiency
	extern const NamedFunc trig_hig;
	extern const NamedFunc eff_higtrig;
}

#endif
