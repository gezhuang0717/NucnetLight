//////////////////////////////////////////////////////////////////////////////
// This file was originally written by Bradley S. Meyer.
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//!
//! \file
//! \brief A header file for user-defined neutrino rate functions.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef USER_NEUTRINO_RATES_H
#define USER_NEUTRINO_RATES_H

#include <boost/unordered_map.hpp>

#include <Libnucnet.h>
#include <WnMatrix.h>
#include <gsl/gsl_const_num.h>
#include <gsl/gsl_const_cgsm.h>

#include "nnt/auxiliary.h"
#include "nnt/math.h"

namespace user
{

//##############################################################################
// Defines
//##############################################################################

#define NU_NUCL                  "neutrino-nucleus"
#define NU_N_CAPTURE             "capture on free neutron"
#define NU_P_CAPTURE             "capture on free proton"
#define S_RHO_RES                "resonant density"
#define S_T_NU                   "Nu T"
#define S_T_NU_E                 "Tnu_e"
#define S_L_NU_E                 "Lnu_e"
#define S_L_NU_E_0               "Lnu_e_0"
#define S_T_NU_E_BAR             "Tnu_e_bar"
#define S_L_NU_E_BAR             "Lnu_e_bar"
#define S_L_NU_E_BAR_0           "Lnu_e_bar_0"
#define S_T_NU_MU                "Tnu_mu"
#define S_L_NU_MU                "Lnu_mu"
#define S_L_NU_MU_0              "Lnu_mu_0"
#define S_T_NU_MU_BAR            "Tnu_mu_bar"
#define S_L_NU_MU_BAR            "Lnu_mu_bar"
#define S_L_NU_MU_BAR_0          "Lnu_mu_bar_0"
#define S_T_NU_TAU               "Tnu_tau"
#define S_L_NU_TAU               "Lnu_tau"
#define S_L_NU_TAU_0             "Lnu_tau_0"
#define S_T_NU_TAU_BAR           "Tnu_tau_bar"
#define S_L_NU_TAU_BAR           "Lnu_tau_bar"
#define S_L_NU_TAU_BAR_0         "Lnu_tau_bar_0"
#define S_LOG10_XSEC             "log_10_xsec"
#define S_PREFACTOR              "rate prefactor"
#define S_DELTA                  "Delta"
#define S_NU_NUCL                "neutrino-nucleus"

//##############################################################################
// Class to encapsulate neutrino arrays.
//##############################################################################

class NeutrinoQuantity
{

  public:
    NeutrinoQuantity( Libnucnet__Reaction * );
    NeutrinoQuantity( const NeutrinoQuantity& );
    ~NeutrinoQuantity();
    double computeValue( double );

  private:
    Libnucnet__Reaction * pReaction;
    std::string sReaction;
    gsl_vector * pTVector;
    gsl_vector * pLog10XSecVector;

};

//##############################################################################
// Prototypes.
//##############################################################################

void
register_my_neutrino_rate_functions( Libnucnet__Reac * );

void
update_my_neutrino_rate_functions_data( nnt::Zone & );

void
set_nu_nucl_hash( Libnucnet__Reac * );

double
nu_nucl_function(
  Libnucnet__Reaction *,
  double,
  nnt::Zone&
);

double
nu_n_capture_function(
  Libnucnet__Reaction *,
  double,
  nnt::Zone&
);

double
nu_p_capture_function(
  Libnucnet__Reaction *,
  double,
  nnt::Zone&
);

int
get_reactant_neutrino(
  Libnucnet__Reaction__Element *, char *
);

double
compute_neutrino_rate(
  nnt::Zone&,
  Libnucnet__Reaction *,
  const char *,
  const char *
);

double
compute_neutrino_cross_section(
  Libnucnet__Reaction *,
  double
);

void
swap_neutrinos( nnt::Zone& );

} // namespace user

#endif // USER_NEUTRINO_RATES_H
