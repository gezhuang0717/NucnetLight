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
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//!
//! \file
//! \brief A header file for routines to write network xml output.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef NNT_WRITE_OUTPUT_XML_H
#define NNT_WRITE_OUTPUT_XML_H

//##############################################################################
// Includes.
//##############################################################################

#include <Libnucnet.h>
#include "nnt/string_defs.h"
#include "nnt/param_defs.h"

namespace nnt
{

  //############################################################################
  // Prototypes.
  //##########################################################################*/

  Libnucnet * create_output( Libnucnet * );

  int add_species_to_zone( Libnucnet__Species *, void * );

  int zone_sort_function( Libnucnet__Zone *, Libnucnet__Zone * );

  void write_xml( Libnucnet *, Libnucnet__Zone * );

  void
  copy_properties(
    const char *,
    const char *,
    const char *,
    const char *,
    Libnucnet__Zone *
  );   

  int insert_species_in_nuc( Libnucnet__Species *, Libnucnet * );

  int insert_reaction_in_reac( Libnucnet__Reaction *, Libnucnet * );

} // namespace nnt
 
#endif /* NNT_WRITE_OUTPUT_XML_H */
