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

	// calculate effect of systematics calculated for each background 
	// in the data control regions on the total bkg. kappa
	extern const NamedFunc wgt_syst_ttx;
	extern const NamedFunc wgt_syst_vjets;
	extern const NamedFunc wgt_syst_qcd;

	// estimate the systematic due to limited knowledge on composition
	extern const NamedFunc wgt_2xhighnb_zjets;

	// analysis trigger and its efficiency
	extern const NamedFunc trig_hig;
	extern const NamedFunc eff_higtrig;
}

#endif
