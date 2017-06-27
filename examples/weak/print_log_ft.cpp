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
//! \brief Example code to print out log10 ft values.
////////////////////////////////////////////////////////////////////////////////

#include <Libnucnet__Reac.h>
#include "nnt/two_d_weak_rates.h"
#include "nnt/auxiliary.h"
#include "nnt/iter.h"

int main( int argc, char * argv[] ) {

  Libnucnet__Reac *p_my_reactions, *p_log10_ft_reactions;

  //============================================================================
  // Check input.
  //============================================================================

  if ( argc != 5 && argc!= 6 ) {
     fprintf(
       stderr, "\nUsage: %s filename t9 rho Ye xpath\n\n", argv[0]
     );
     fprintf(
       stderr, "  filename = input xml rate data filename\n\n"
     );
     fprintf(
       stderr, "  t9 = input temperature\n\n"
     );
     fprintf(
       stderr, "  rho = input density (g/cc)\n\n"
     );
     fprintf(
       stderr, "  Ye = input electron-to-baryon ratio (g/cc)\n\n"
     );
     fprintf(
       stderr, "  xpath = xpath expression (optional)\n\n"
     );

     return EXIT_FAILURE;
  }

  //============================================================================
  // Read file and exit if not present.
  //============================================================================

  if ( argc == 4 )
    p_my_reactions = Libnucnet__Reac__new_from_xml( argv[1], NULL );
  else
    p_my_reactions = Libnucnet__Reac__new_from_xml( argv[1], argv[5] );

  //============================================================================
  // Get log ft reactions.
  //============================================================================

  p_log10_ft_reactions =
    Libnucnet__Reac__extractSubset(
      p_my_reactions,
      "[user_rate/@key = 'two-d weak rates log10 ft']"
    ); 

  //============================================================================
  // Done with p_my_reactions.
  //============================================================================

  Libnucnet__Reac__free( p_my_reactions );

  //============================================================================
  // Print log ft.
  //============================================================================

  fprintf( stdout, "\n\t\t\t\tReaction\t\t\tlog10 ft\n" );
  fprintf(
    stdout,
    "=========================================================\t============\n"
  );

  Libnucnet__Reac__setReactionCompareFunction(
    p_log10_ft_reactions,
    (Libnucnet__Reaction__compare_function) nnt::compare_reactions_by_string
  );

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list( p_log10_ft_reactions );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    nnt::TwoDWeakQuantity
      my_log10_ft( reaction.getNucnetReaction(), nnt::s_LOG10_FT );

    fprintf(
      stdout,
      "%-57s\t%12.6e\n",
      Libnucnet__Reaction__getString( reaction.getNucnetReaction() ),
      my_log10_ft.computeValue(
        atof( argv[2] ),
        atof( argv[3] ) * atof( argv[4] )
      ).first
    );

  }

  fprintf( stdout, "\n" );

  //============================================================================
  // Print number of reactions.
  //============================================================================

  fprintf(
    stdout,
    "Number of reactions = %lu\n\n",
    (unsigned long)
       Libnucnet__Reac__getNumberOfReactions( p_log10_ft_reactions )
  );

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__Reac__free( p_log10_ft_reactions );

  return EXIT_SUCCESS;

}

