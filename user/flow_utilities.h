////////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
//!
//! \file
//! \brief A header file for routines related to reaction flows.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef USER_FLOW_UTILITIES_H
#define USER_FLOW_UTILITIES_H

#include <iostream>
#include <utility>

#include <Libnucnet.h>
#include <gsl/gsl_vector.h>

#include "nnt/auxiliary.h"
#include "nnt/string_defs.h"
#include "nnt/iter.h"
#include "nnt/weak_detailed_balance.h"

#include "user/rate_modifiers.h"
#include "user/user_rate_functions.h"
#include "user/coul_corr.h"

/**
 * @brief A namespace for user-defined rate functions.
 */
namespace user
{

std::pair<double,double>
compute_flows_for_reaction( nnt::Zone&, Libnucnet__Reaction * );

std::pair<double,double>
compute_modified_flows_for_reaction(
  nnt::Zone&,
  Libnucnet__Reaction *
);

gsl_vector *
compute_flow_vector( nnt::Zone&, Libnucnet__NetView * );

gsl_vector *
compute_modified_flow_vector( nnt::Zone&, Libnucnet__NetView * );

gsl_vector *
compute_forward_flow_vector( nnt::Zone&, Libnucnet__NetView * );

gsl_vector *
compute_reverse_flow_vector( nnt::Zone&, Libnucnet__NetView * );

gsl_vector *
compute_modfied_flow_vector( nnt::Zone&, Libnucnet__NetView * );

double
compute_total_flow(
  nnt::Zone&,
  const char *,
  Libnucnet__NetView *
);

double
compute_entropy_change_rate(
  nnt::Zone&,
  Libnucnet__NetView *
);

double
compute_energy_generation_rate_per_nucleon(
  nnt::Zone&,
  Libnucnet__NetView *
);

void
update_flow_currents( nnt::Zone&, nnt::Zone& );

} // namespace user

#endif // USER_FLOW_UTILITIES_H

