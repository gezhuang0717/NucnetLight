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
//! \brief A header file defining functions to remove duplicate reactions.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef USER_REMOVE_DUPLICATES_H
#define USER_REMOVE_DUPLICATES_H

/*##############################################################################
// Includes.
//############################################################################*/

#include <Libnucnet.h>

namespace user
{

/*##############################################################################
// Define. Change "no" to "yes" to print out information about removed
// reactions.
//############################################################################*/

#define PRINTOUT "yes"

/*##############################################################################
// Prototypes.
//############################################################################*/

void
remove_duplicate_reactions( Libnucnet__Net * );

int
remove_duplicate( Libnucnet__Reaction *, Libnucnet__Net * );

int
get_number_reactants( Libnucnet__Reaction__Element *, int * );

} // namespace user

#endif // USER_REMOVE_DUPLICATES_H
