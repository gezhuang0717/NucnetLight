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
//! \file
//! \brief A file containing routines to write network xml output.
//!
////////////////////////////////////////////////////////////////////////////////

#include "write_output_xml.h"

namespace nnt
{

  //############################################################################
  // create_output()
  //############################################################################

  Libnucnet *
  create_output(
    Libnucnet *self
  )
  {

    Libnucnet *p_new_nucnet;

    p_new_nucnet =  Libnucnet__new();

    Libnucnet__Nuc__iterateSpecies(
      Libnucnet__Net__getNuc( Libnucnet__getNet( self ) ),
      (Libnucnet__Species__iterateFunction) insert_species_in_nuc,
      p_new_nucnet
    );

    Libnucnet__Reac__iterateReactions(
      Libnucnet__Net__getReac( Libnucnet__getNet( self ) ),
      (Libnucnet__Reaction__iterateFunction) insert_reaction_in_reac,
      p_new_nucnet
    );

    return p_new_nucnet;

  }

  //############################################################################
  // write_xml()
  //############################################################################

  void
  write_xml(
    Libnucnet *p_zone_data,
    Libnucnet__Zone *p_zone
  )
  {

    Libnucnet__Zone * p_new_zone;
    gsl_vector * p_vector;
    char s_zone[i_BUF_SIZE];

    sprintf(
      s_zone,
      "%lu",
      (unsigned long) Libnucnet__getNumberOfZones( p_zone_data ) + 1
    );

    p_new_zone =
      Libnucnet__Zone__new(
        Libnucnet__getNet( p_zone_data ),
        s_zone,
        NULL,
        NULL
      );

    Libnucnet__Zone__iterateOptionalProperties(
      p_zone,
      NULL,
      NULL,
      NULL,
      (Libnucnet__Zone__optional_property_iterate_function)
         copy_properties,
      p_new_zone
    );

    p_vector = Libnucnet__Zone__getAbundances( p_zone );
    Libnucnet__Zone__updateAbundances( p_new_zone, p_vector );
    gsl_vector_free( p_vector );

    p_vector = Libnucnet__Zone__getAbundanceChanges( p_zone );
    Libnucnet__Zone__updateAbundanceChanges( p_new_zone, p_vector );
    gsl_vector_free( p_vector );

    if( !Libnucnet__addZone( p_zone_data, p_new_zone ) )
    {
      fprintf( stderr, "Couldn't add zone!\n" );
      exit( EXIT_FAILURE );
    } 

  }

  //############################################################################
  // copy_properties
  //############################################################################

  void
  copy_properties(
    const char * s_name,
    const char * s_tag1,
    const char * s_tag2,
    const char * s_value,
    Libnucnet__Zone * p_zone
  )
  {

    Libnucnet__Zone__updateProperty(
      p_zone,
      s_name,
      s_tag1,
      s_tag2,
      s_value
    );

  } 

  //############################################################################
  // zone_sort_function()
  //############################################################################

  int
  zone_sort_function( Libnucnet__Zone *p_zone1, Libnucnet__Zone *p_zone2 )
  {

    if(
       atoi( Libnucnet__Zone__getLabel( p_zone1, 1 ) ) <
       atoi( Libnucnet__Zone__getLabel( p_zone2, 1 ) )
    )
       return -1;
    else if(
       atoi( Libnucnet__Zone__getLabel( p_zone1, 1 ) ) >
       atoi( Libnucnet__Zone__getLabel( p_zone2, 1 ) )
    )
       return 1;
    else
       return 0;

  }

  //############################################################################
  // insert_species_in_nuc().
  //############################################################################

  int
  insert_species_in_nuc(
    Libnucnet__Species *p_species,
    Libnucnet *p_nucnet
  )
  {

    Libnucnet__Nuc__addSpecies(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_nucnet ) ),
      Libnucnet__Species__copy( p_species )
    );

    return 1;

  }

  //############################################################################
  // insert_reaction_in_reac().
  //############################################################################

  int
  insert_reaction_in_reac(
    Libnucnet__Reaction *p_reaction,
    Libnucnet *p_nucnet
  )
  {

    if(
       Libnucnet__Net__isValidReaction(
	 Libnucnet__getNet( p_nucnet ),
	 p_reaction
       )
    )
      Libnucnet__Reac__addReaction(
	Libnucnet__Net__getReac( Libnucnet__getNet( p_nucnet ) ),
	Libnucnet__Reaction__copy( p_reaction )
      );

    return 1;

  }

}
