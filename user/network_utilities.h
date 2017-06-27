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
//! \brief A header file to define useful utilities for network routines.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef NETWORK_UTILITIES_H
#define NETWORK_UTILITIES_H

#include <assert.h>

#include "nnt/auxiliary.h"
#include "nnt/math.h"
#include "nnt/string_defs.h"
#include "nnt/wrappers.h"

/**
 * @brief A namespace for user-defined rate functions.
 */
namespace user
{

//##############################################################################
// Numerical defines.
//##############################################################################

#define  D_DT_MIN  1.e-20
#define  D_EPS     0.05
#define  D_GAMMA   5. / 3.

//##############################################################################
// Prototypes.
//##############################################################################

int set_zone_post_shock_conditions( nnt::Zone, std::string );

double compute_luminosity( nnt::Zone );

double compute_photon_entropy_loss_rate( nnt::Zone );

double compute_v_Thermal( Libnucnet__Species *, double );

void update_exposures( nnt::Zone& );

void
update_exposures_callback(
  const char *,
  const char *,
  const char *,
  const char *,
  nnt::Zone&
);

double compute_Ysum( nnt::Zone&, nnt::species_list_t& );

void
update_t9_rho_in_zone_by_interpolation(
  nnt::Zone&,
  std::string,
  const std::vector<double>&,
  const std::vector<double>&,
  const std::vector<double>&
);

void
copy_zone_abundances_as_properties(
  nnt::Zone&,
  nnt::Zone&,
  const char *
);

double
compute_cluster_abundance_moment(
  nnt::Zone&,
  const char *,
  std::string,
  int
);

void
set_specific_species(
  WnMatrix *,
  gsl_vector *,
  void *
);

} // namespace user

#endif // NETWORK_UTILITIES_H
