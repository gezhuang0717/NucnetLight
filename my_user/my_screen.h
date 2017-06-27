//////////////////////////////////////////////////////////////////////////////
// This file was originally written by Bradley S. Meyer and George C.
// Jordan, IV.
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
//! \brief A header file defining screening functions.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef USER_SCREEN_H
#define USER_SCREEN_H

#include <Libnucnet.h>

#include "nnt/wrappers.h"

namespace my_user
{

/*##############################################################################
// Defines.
//############################################################################*/

#define D_THETA_E   0.0   /* setting to 1.0 for non-degenerate matter.
                             Set to 0.0 for the completely degenerate case.
                          */

/*##############################################################################
// User-supplied data structure.
//############################################################################*/

typedef struct screening_data {
  double dYe2;
} screening_data;

/*##############################################################################
// Prototypes for user-supplied functions.
//############################################################################*/

screening_data
get_screening_data( nnt::Zone& zone );

double
my_screening_function(
  double, double, double, int, int, int, int, screening_data *
);

double calculate_gamma_e( double, double, double );

double calculate_gamma_effective( int, int, double );

double
strong_screening_factor(
  int, int, int, int, double, double, double
);

double
weak_screening_factor(
  int, int, double, double, double, double, double
);

double intermediate_screening_factor( double, double );

} // namespace my_user

#endif // USER_SCREEN_H
