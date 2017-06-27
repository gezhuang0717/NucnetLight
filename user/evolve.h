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
//! \file
//! \brief A header file for the evolution routines.
////////////////////////////////////////////////////////////////////////////////

#ifndef EVOLVE_H
#define EVOLVE_H

//##############################################################################
// Includes. 
//##############################################################################

#include <Libnucnet.h>
#include <Libnuceq.h>
#include <Libstatmech.h>

#include <boost/math/special_functions/erf.hpp>

#include "nnt/auxiliary.h"
#include "nnt/weak_detailed_balance.h"

#include "user/remove_duplicate.h"
#include "user/screen.h"
#include "user/thermo.h"
#include "user/coul_corr.h"
#include "user/network_limiter.h"
#include "user/matrix_solver.h"
#include "user/weak_utilities.h"
#include "user/rate_modifiers.h"


/**
 * @brief A namespace for user-defined rate functions.
 */
namespace user
{

//##############################################################################
// Define some parameters.
//##############################################################################

#define D_MIN          1.e-04  // Newton-Raphson convergence criterion
#define D_X_EPS        1.e-08  // Xsum check for safe evolve
#define D_Y_MIN        1.e-10  // Smallest y for convergence
#define D_Y_MIN_PRINT  0.e-30  // Smallest y to print out
#define I_ITMAX        30      // Maximum number of Newton-Raphson iterations

//##############################################################################
// Enumeration.
//##############################################################################

enum solvers { ARROW, GSL };

//##############################################################################
// Prototypes.
//##############################################################################

int
evolve( nnt::Zone& );

void safe_evolve( nnt::Zone&, double, double );

bool is_nonneg_abunds( nnt::Zone& );

int
evolve_zone( nnt::Zone&, double );

double
compute_neutrino_entropy_loss_rate( nnt::Zone );

void
evolve_nse_plus_weak_rates( nnt::Zone );

double
nse_plus_weak_rates_function( double, void * );

void
set_zone_for_evolution( nnt::Zone& );

WnMatrix *
get_evolution_matrix( nnt::Zone& );

std::pair< WnMatrix *, gsl_vector * >
get_evolution_matrix_and_vector( nnt::Zone& );

std::pair<double,double>
check_matrix_solution(
  nnt::Zone&,
  gsl_vector *
);

double my_log10_t9_from_entropy_root( double, nnt::Zone * );

double my_density_from_entropy_root( double, nnt::Zone * );

std::pair<double,double>
my_log10_t9_guess( nnt::Zone * );

} // namespace user

#endif // EVOLVE_H
