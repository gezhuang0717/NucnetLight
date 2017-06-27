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
//! \brief A header file for the zone matrix solver routines.
////////////////////////////////////////////////////////////////////////////////

#ifndef MATRIX_SOLVER_H
#define MATRIX_SOLVER_H

//##############################################################################
// Includes. 
//##############################################################################

#include <iostream>

#include <boost/unordered_map.hpp>
#include <boost/lexical_cast.hpp>

#include "nnt/wrappers.h"
#include "nnt/string_defs.h"

#ifdef SPARSKIT2
#include <WnSparseSolve.h>
#include "ilu_solvers.h"
#endif

/**
 * @brief A namespace for user-defined functions.
 */
namespace user
{

gsl_vector *
solve_matrix_for_zone( nnt::Zone&, WnMatrix *, gsl_vector * );

#ifdef SPARSKIT2
gsl_vector *
solve_sparse_matrix_with_ilu_preconditioner(
  WnMatrix *,
  gsl_vector *,
  boost::unordered_map<std::string, std::string>&
);

gsl_vector *
phi__solve__parallel(
  WnSparseSolve__Phi *,
  std::vector<WnMatrix *>&,
  const gsl_vector *,
  const gsl_vector *,
  double
);
#endif

} // namespace user

#endif // MATRIX_SOlVER_H
