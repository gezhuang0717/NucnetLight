/*//////////////////////////////////////////////////////////////////////////////
// <file type="public">
//   <license>
//     This file was originally written by Bradley S. Meyer.
//
//     This is free software; you can redistribute it and/or modify it
//     under the terms of the GNU General Public License as published by
//     the Free Software Foundation; either version 2 of the License, or
//     (at your option) any later version.
//
//     This software is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     Please see the src/README.txt file in this distribution for more
//     information.

//   </license>
//   <description>
//     <abstract>
//       Header file for the remove_duplicate code.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

/*##############################################################################
// Includes.
//############################################################################*/

#include <Libnucnet.h>

/*##############################################################################
// Use extern "C" for C++ compilers.
//############################################################################*/

#ifdef __cplusplus
extern "C"
{
#endif

/*##############################################################################
// Define. Change "no" to "yes" to print out information about removed
// reactions.
//############################################################################*/

#define PRINTOUT "no"

/*##############################################################################
// Prototype.
//############################################################################*/

int
remove_duplicate( Libnucnet__Reaction *, Libnucnet__Net * );
#ifdef __cplusplus
} /* extern "C" */
#endif
