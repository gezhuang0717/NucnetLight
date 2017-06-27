// Copyright (c) 2011-2015 Clemson University.
//
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
// You should have received a copy of the GNU General Public License
// along with this software; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
// USA
//
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//!
//! \file
//! \brief A header file to define NucNet Tools strings.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef NNT_STRING_DEFS_H
#define NNT_STRING_DEFS_H

namespace nnt
{

  const char s_ANTI_NEUTRINO_E[]              = "anti-neutrino_e";
  const char s_ARROW[]                        = "Arrow";
  const char s_ARROW_WIDTH[]                  = "Arrow width";
  const char s_AVERAGE_ENERGY_NU_E[]          = "av nu_e energy";
  const char s_AVERAGE_ENERGY_NUBAR_E[]       = "av nubar_e energy";
  const char s_BARYON[]                       = "baryon";
  const char s_BARYON_ENTROPY[]               = "s_b";
  const char s_BETA_MINUS_XPATH[]             = \
     "[product = 'electron' and product = 'anti-neutrino_e']";
  const char s_BETA_PLUS_XPATH[]              = \
     "[product = 'positron' and product = 'neutrino_e']"; 
  const char s_CLUSTER[]                      = "cluster";
  const char s_CLUSTER_XPATH[]                = "cluster_xpath";
  const char s_CLUSTER_CONSTRAINT[]           = "cluster_constraint";
  const char s_CLUSTER_CHEMICAL_POTENTIAL[]   = "muh_kT";
  const char s_DMUEKTDT[]                     = "d_dmuekT_dT";
  const char s_DPDT[]                         = "dPdT";
  const char s_DTIME[]                        = "dt";
  const char s_DTMIN[]                        = "dtmin";
  const char s_DTMAX[]                        = "dtmax";
  const char s_DT_SAFE[]                      = "dt safe";
  const char s_ELECTRON[]                     = "electron";
  const char s_ELECTRON_CAPTURE_XPATH[]       = \
    "[reactant = 'electron' and product = 'neutrino_e']"; 
  const char s_ELECTRON_ENTROPY[]             = "s_e";
  const char s_ENTROPY_PER_NUCLEON[]          = "entropy per nucleon";
  const char s_EQUIL_ONLY[]                   = "only equilibrium";
  const char s_EXPOSURE[]                     = "exposure";
  const char s_FINAL_ABUNDANCE[]              = "final abundance";
  const char s_FLOW_CURRENT[]                 = "flow current";
  const char s_FORWARD_FLOW[]                 = "foward";
  const char s_GSL[]                          = "Gsl";
  const char s_HI_T9_EQUIL[]                  = "high t9 for equilibrium";
  const char s_ILU_DELTA[]                    = "ilu delta";
  const char s_ILU_DROP_TOL[]                 = "ilu drop tolerance";
  const char s_INITIAL_ABUNDANCE[]            = "initial abundance";
  const char s_INTERNAL_ENERGY_DENSITY[]      = "internal energy density";
  const char s_ITER_SOLVER[]                  = "iterative solver method";
  const char \
    s_ITER_SOLVER_CONVERGENCE_METHOD[]        = \
        "iterative solver convergence method";
  const char s_ITER_SOLVER_DEBUG[]            = "iterative solver debug";
  const char s_ITER_SOLVER_MAX_ITERATIONS[]   = \
        "iterative solver maximum iterations";
  const char s_ITER_SOLVER_REL_TOL[]   = \
        "iterative solver relative tolerance";
  const char s_ITER_SOLVER_ABS_TOL[]   = \
        "iterative solver absolute tolerance";
  const char s_ITER_SOLVER_T9[]               = "t9 for iterative solver";
  const char s_LAB_RATE[]                     = "lab rate";
  const char s_LAB_RATE_T9_CUTOFF[]           = "lab rate t9 cutoff";
  const char s_LAB_RATE_T9_CUTOFF_FACTOR[]    = \
    "lab rate t9 cutoff factor";
  const char s_LARGE_NEG_ABUND_THRESHOLD[]    = \
    "large negative abundances threshold";
  const char s_LOG10_FT[]                     = "log10_ft";
  const char s_LOG10_RATE[]                   = "log10_rate";
  const char s_LOG10_RHOE[]                   = "log10_rhoe";
  const char s_LOW_T_NETWORK_T9[]             = "low T network T9";
  const char s_LOW_T_NETWORK_REAC_XPATH[]     = "low T network reaction xpath";
  const char s_MACH[]                         = "mach number";
  const char s_MODIFIED_NET_FLOW[]            = "modified net flow";
  const char s_MUEKT[]                        = "muekT";
  const char s_MUPKT[]                        = "mupkT";
  const char s_MUNKT[]                        = "munkT";
  const char s_MU_NUE_KT[]                    = "munuekT";
  const char s_MY_ENTROPY_DENSITY[]           = "my_total_entropy_density";
  const char s_MY_FLOW_MIN[]                  = "my minimum net flow";
  const char s_NET_FLOW[]                     = "net";
  const char s_NEUTRINO_E[]                   = "neutrino_e";
  const char s_NEUTRINO_REACTION_XPATH[]      = \
    "[reactant[contains(.,'neutrino')]]";
  const char s_NEWTON_RAPHSON_ABUNDANCE[]     = \
    "Newton-Raphson abundance minimum";
  const char s_NEWTON_RAPHSON_CONVERGE[]      = \
    "Newton-Raphson convergence minimum";
  const char s_NUC_ENERGY[]                   = "nuc energy";
  const char s_NUMBER_CLUSTERS[]              = "clusters";
  const char s_OLD_T9[]                       = "old t9";
  const char s_OLD_RHO[]                      = "old rho";
  const char s_OLD_ZONE[]                     = "old_zone";
  const char s_PHOTON[]                       = "photon";
  const char s_PHOTON_ENTROPY[]               = "s_p";
  const char s_POSITRON[]                     = "positron";
  const char s_POSITRON_CAPTURE_XPATH[]       = \
    "[reactant = 'positron' and product = 'anti-neutrino_e']"; 
  const char s_PRESSURE[]                     = "pressure";
  const char s_RADIOACTIVE_MASS_FRACTION[]    = "radioactive mass";
  const char s_RADIUS_0[]                     = "radius_0";
  const char s_RADIUS[]                       = "radius";
  const char s_RATE_MODIFICATION_VIEW[]       = "rate modification view";
  const char s_RATE_MODIFICATION_FACTOR[]     = "factor"; 
  const char s_RATE_MODIFICATION_NUC_XPATH[]  = "nuclide xpath";
  const char s_RATE_MODIFICATION_REAC_XPATH[] = "reaction xpath";
  const char s_RATE_THRESHOLD[]               = "rate threshold";
  const char s_REVERSE_FLOW[]                 = "reverse";
  const char s_RHO[]                          = "rho";
  const char s_RHO_0[]                        = "rho_0";
  const char s_RHO1[]                         = "cell density";
  const char s_SDOT[]                         = "sdot";
  const char s_SOLVER[]                       = "solver";
  const char s_SHOCK_SPEED[]                  = "shock speed";
  const char s_SMALL_ABUNDANCES_THRESHOLD[]   = "small abundances threshold";
  const char s_SMALL_RATES_THRESHOLD[]        = "small rates threshold";
  const char s_SOUND_SPEED[]                  = "sound speed";
  const char s_SPECIFIC_ABUNDANCE[]           = "specific abundance";
  const char s_SPECIFIC_HEAT_PER_NUCLEON[]    = "cv";
  const char s_SPECIFIC_HEAT_PER_NUCLEON_AT_CONSTANT_PRESSURE[] \
                                              = "cp";
  const char s_SPECIES_REMOVAL_NUC_XPATH[]   = "species removal nuclide xpath";
  const char s_SPECIES_REMOVAL_REAC_XPATH[]   =\
    "species removal reaction xpath";
  const char s_SPECIFIC_SPECIES[]             = "specific species";
  const char s_STEPS[]                        = "steps";
  const char s_T1[]                           = "cell temperature";
  const char s_TEND[]                         = "tend";
  const char s_TAU[]                          = "tau";
  const char s_TAU_LUM_NEUTRINO[]             = "tau for neutrino luminosity";
  const char s_THERMO_NUC_VIEW[]              = "thermo nuc view";
  const char s_TIME[]                         = "time";
  const char s_TOTAL[]                        = "total";
  const char s_TWO_D_WEAK_RATES[]             = "two-d weak rates";
  const char s_TWO_D_WEAK_RATES_LOG10_FT[]    = "two-d weak rates log10 ft";
  const char s_TWO_D_WEAK_XPATH []            = \
    "[user_rate/@key = 'two-d weak rates log10 ft' or \
      user_rate/@key = 'two-d weak rates']";
  const char s_T9[]                           = "t9";
  const char s_T9_0[]                         = "t9_0";
  const char s_USE_HI_T_EQUIL[]               = \
        "use high temperature equilibrium";
  const char s_USE_SCREENING[]                = "use screening";
  const char s_USE_APPROXIMATE_WEAK_RATES[]   = "use approximate weak rates";
  const char s_USE_WEAK_DETAILED_BALANCE[]    = "use weak detailed balance";
  const char s_WEAK_VIEW_FOR_LAB_RATE_TRANSITION[] = \
    "weak view for lab rate transition";
  const char s_WEAK_XPATH[]                   = \
    "[reactant = 'electron' or product = 'electron' or \
      reactant = 'positron' or product = 'positron']";
  const char s_YC[]                           = "Yc";
  const char s_YCDOT[]                        = "Ycdot";
  const char s_YE[]                           = "Ye";
  const char s_YEDOT[]                        = "Yedot";
  const char s_YEDOT_BM[]                     = "Yedot_bm";
  const char s_YEDOT_EC[]                     = "Yedot_ec";
  const char s_YEDOT_PC[]                     = "Yedot_pc";
  const char s_YEDOT_BP[]                     = "Yedot_bp";
  const char s_ZONE_LABEL[]                   = "zone label";
  const char s_ZONE_MASS[]                    = "zone mass";
  const char s_ZONE_MASS_CHANGE[]             = "zone mass change";

} // namespace nnt

#endif /* NNT_STRING_DEFS_H */
