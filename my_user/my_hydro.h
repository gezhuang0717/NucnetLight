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
//! \brief A header file for the thermodynamics evolution file.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <iostream>
#include <fstream>

#include "nnt/auxiliary.h"
#include "nnt/string_defs.h"
#include "nnt/math.h"

#include "user/network_utilities.h"

/**
 * @brief A namespace for user-defined rate functions.
 */
namespace my_user
{

//##############################################################################
// Defines.
//##############################################################################

#define S_OUTPUT_FILE   "output file"

//##############################################################################
// Validation.  "no" = no validation, "yes" = validation.
//##############################################################################

#define VALIDATE       "no"

//##############################################################################
// Prototypes.
//##############################################################################

Libnucnet * get_nucnet( int, char ** );

void
get_trajectory_data( char * );

void
initialize_zone( nnt::Zone&, char ** );

void update_zone_properties( nnt::Zone& );

int set_zone( Libnucnet *, nnt::Zone&, char ** );

} // namespace my_user
