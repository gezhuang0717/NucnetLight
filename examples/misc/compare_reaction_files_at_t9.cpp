////////////////////////////////////////////////////////////////////////////////
// This file was originally written by Tianhong Yu and Bradley S. Meyer.
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
//! \brief Example code to compare reaction rates at a given temperature
//          from two different files.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <gsl/gsl_math.h>
#include <Libnucnet__Reac.h>
#include "nnt/auxiliary.h"

//##############################################################################
// Defines.
//##############################################################################

#define I_SHOW_IF_NOT_PRESENT 0 // 1 = show reactions in file 1 not in file 2.

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  Libnucnet__Reac * p_reac_1 = NULL, * p_reac_2;
  Libnucnet__Reaction * p_reaction;
  double d_rate_1, d_rate_2;

  //============================================================================
  // Check input.
  //============================================================================

  if ( argc < 5 || argc > 6 ) {
    fprintf(
      stderr,
      "\nUsage: %s xml_file1 xml_file2 t9 xpath_reac\n\n",
      argv[0]
    );
    fprintf(
      stderr, "  xml_file1 = first reaction data xml filename\n\n"
    );
    fprintf(
      stderr, "  xml_file2 = second reaction data xml filename\n\n"
    );
    fprintf(
      stderr, "  t9 = temperature at which to compare rates\n\n"
    );
    fprintf(
      stderr, "  factor = factor by which to compare rates\n\n"
    );
    fprintf(
      stderr,
      "  xpath_reac = reaction xpath expression (optional)\n\n"
    );
    return EXIT_FAILURE;
  }

  //============================================================================
  // Read and store input.
  //============================================================================

  if( argc == 5 )
  {
    p_reac_1 =
      Libnucnet__Reac__new_from_xml(
        argv[1],
        NULL
      );  
    p_reac_2 =
      Libnucnet__Reac__new_from_xml(
        argv[2],
        NULL
      );  
  }
  else if( argc == 6 )
  {
    p_reac_1 =
      Libnucnet__Reac__new_from_xml(
        argv[1],
        argv[5]
      );  
    p_reac_2 =
      Libnucnet__Reac__new_from_xml(
        argv[2],
        argv[5]
      );  
  }

  //============================================================================
  // Get reaction list.
  //============================================================================

  Libnucnet__Reac__setReactionCompareFunction(
    p_reac_1,
    (Libnucnet__Reaction__compare_function) nnt::compare_reactions_by_string
  );

  nnt::reaction_list_t reaction_list = nnt::make_reaction_list( p_reac_1 );

  //============================================================================
  // Print header.
  //============================================================================

  fprintf(
    stdout,
    "\n\t\t\tReaction\t\t\t\tFile 1          File 2\n"
  );
  fprintf(
    stdout,
    "=======================================================         ============    ============\n"
  );

  //============================================================================
  // Loop.
  //============================================================================

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    p_reaction =
      Libnucnet__Reac__getReactionByString(
        p_reac_2,
        Libnucnet__Reaction__getString( reaction.getNucnetReaction() )
      );

    if( p_reaction )
    {

      d_rate_1 =
        Libnucnet__Reaction__computeRate(
          reaction.getNucnetReaction(),
          atof( argv[3] ),
          NULL
        );

      d_rate_2 =
        Libnucnet__Reaction__computeRate(
          p_reaction,
          atof( argv[3] ),
          NULL
        );

      if( gsl_fcmp( d_rate_1, d_rate_2, atof( argv[4] ) ) )
      {
        fprintf(
          stdout,
          "%-57s\t%12.6e\t%12.6e\n",
          Libnucnet__Reaction__getString( reaction.getNucnetReaction() ),
          d_rate_1,
          d_rate_2
        );
      }
     
    }
    else
    {

      if( I_SHOW_IF_NOT_PRESENT )
      {
        fprintf(
          stdout,
          "%-57s\t%12.6e\t%s\n",
          Libnucnet__Reaction__getString( reaction.getNucnetReaction() ),
          d_rate_1,
          "---------"
        );
      }

    }

  }

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__Reac__free( p_reac_1 );
  Libnucnet__Reac__free( p_reac_2 );

  return EXIT_SUCCESS;

}

