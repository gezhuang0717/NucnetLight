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
//! \brief A header for user-defined Coulomb correction functions.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef USER_COUL_CORR_H
#define USER_COUL_CORR_H

#include <Libnucnet.h>

namespace user
{

//##############################################################################
// User-supplied data structure.
//##############################################################################

typedef struct coul_corr_data{
  double dA1;
  double dA2;
} coul_corr_data;


//##############################################################################
// Prototypes for user-supplied function.
//##############################################################################

coul_corr_data
get_coulomb_corr_data( void );

double
my_coulomb_correction(
  Libnucnet__Species *, double, double, double, coul_corr_data *
);

} // namespace user

#endif // USER_COUL_CORR_H

