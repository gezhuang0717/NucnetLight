////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011 Clemson University.
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
//! \file
//! \brief A header file to define NucNet Tools thermodynamics routines.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef THERMO_HPP
#define THERMO_HPP

#include <iostream>

#include <Libnucnet.h>
#include <Libstatmech.h>

#include "nnt/string_defs.h"
#include "nnt/param_defs.h"
#include "nnt/iter.h"
#include "nnt/auxiliary.h"

namespace user
{

nnt::species_list_t
get_thermo_species_list( nnt::Zone& );

double
compute_log10_t9_entropy_root_with_equilibrium(
  double,
  nnt::Zone *
);

double
compute_thermo_quantity(
  nnt::Zone&,
  std::string,
  std::string
);

double
compute_sound_speed( nnt::Zone& );

double
compute_log10_t9_entropy_root(
  double,
  nnt::Zone *
);

double
compute_log10_density_entropy_root_with_equilibrium(
  double,
  nnt::Zone *
);

double
compute_dlnG_dT( Libnucnet__Species *, double );

double
compute_lnG( double, void * );

double
dP_drho( double, void * );

std::pair<double, double>
sound_speed_t9_guess_function( nnt::Zone * );

double
dS_dT_at_constant_pressure(
  double,
  void *
);

double
compute_log10_rho_pressure_root(
  double,
  void *
);

std::pair<double,double>
cp_rho_guess_function( nnt::Zone * );

} // namespace user

#endif // THERMO_HPP
