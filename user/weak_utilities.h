////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2012 Clemson University.
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
//////////////////////////////////////////////////////////////////////////////*/

////////////////////////////////////////////////////////////////////////////////
//!
//! \file weak_utilities.h
//! \brief A header file to define useful utilities for weak reactions.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef WEAK_UTILITIES_H
#define WEAK_UTILITIES_H

#include <boost/unordered_map.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/math/special_functions/erf.hpp>

#include "nnt/auxiliary.h"
#include "nnt/string_defs.h"
#include "nnt/wrappers.h"
#include "nnt/two_d_weak_rates.h"

#include "user/user_rate_functions.h"
#include "user/flow_utilities.h"

namespace user
{

//##############################################################################
// Defines.
//##############################################################################

#define D_EPS_ABS               0.
#define D_EPS_REL               1.e-12

#define S_LAMBDA_3              "lambda 3"
#define S_DYEDOT_DYE            "dyedot_dye"
#define S_DYEDOT_DYE_NUCLEON    "dyedot_dye_nucleon"

//##############################################################################
// Prototypes.
//##############################################################################

double
my_yedot(
  double,
  void *
);

boost::tuple<double, double, double, double, double>
compute_all_yedot(
  double,
  void *
);

double
compute_dyedot_dye(
  nnt::Zone&,
  double
);

void
set_weak_views_in_zones( Libnucnet * );

double
my_yedot_root( double, void * );

double
compute_lambda( nnt::Zone );

double compute_approximate_neutrino_entropy_loss_rate( nnt::Zone & );

double compute_neutrino_entropy_loss_rate( nnt::Zone );

double compute_neutrino_energy_loss_rate( nnt::Zone );

double 
compute_reaction_neutrino_energy_loss_rate(
  Libnucnet__Reaction *,
  Libnucnet__Net *,
  double, 
  double,
  double,
  double,
  double
);

double 
compute_simple_neutrino_entropy_loss_rate(
  nnt::Zone,
  double
);

int
set_neutrino_type(
  Libnucnet__Reaction__Element *,
  char *
);

double 
compute_weak_rate_from_log10_ft(
  Libnucnet__Reaction *,
  Libnucnet__Net *,
  double,
  double,
  double,
  double,
  double
);

void
set_two_d_weak_rates_hashes( Libnucnet__Reac * );

void
set_two_d_weak_energy_loss_hash( Libnucnet__Reac * );

double
my_log10_ft_rate_function(
  Libnucnet__Reaction *,
  double,
  void *
);

double
compute_two_d_weak_rate(
  Libnucnet__Reaction *,
  double,
  double *
);

double
my_approximate_weak_rate_function(
  Libnucnet__Reaction *,
  double,
  void *
);

void
correct_for_weak_lab_rate(
  Libnucnet__Reaction *,
  double,
  double &
);

} // namespace user

#endif // WEAK_UTILITIES_H
