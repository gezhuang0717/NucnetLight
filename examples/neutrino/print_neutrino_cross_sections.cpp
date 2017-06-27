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
//! \brief Example code to print neutrino-nucleus cross sections.
////////////////////////////////////////////////////////////////////////////////

#include <Libnucnet__Reac.h>
#include "nnt/string_defs.h"
#include "nnt/auxiliary.h"
#include "nnt/iter.h"

#include "user/neutrino_rate_functions.h"

int main( int argc, char * argv[] ) {

  Libnucnet__Reac *p_my_reactions;
  Libnucnet__ReacView * p_view;
  std::string s_xpath;

  //============================================================================
  // Check input.
  //============================================================================

  if ( argc != 3 && argc!= 4 ) {
     fprintf(
       stderr, "\nUsage: %s neutrino_xml xpath T_nu xpath\n\n", argv[0]
     );
     fprintf(
       stderr, "  neutrino_xml = input neutrino data xml filename\n\n"
     );
     fprintf(
       stderr, "  T_nu = input neutrino temperature\n\n"
     );
     fprintf(
       stderr, "  xpath = reaction xpath expression (optional)\n\n"
     );
     return EXIT_FAILURE;
  }

  //============================================================================
  // Read file and exit if not present.
  //============================================================================

  if ( argc == 3 )
    p_my_reactions =
      Libnucnet__Reac__new_from_xml(
        argv[1],
        NULL
      );
  else
    p_my_reactions =
      Libnucnet__Reac__new_from_xml(
        argv[1],
        argv[3]
    );

  if( !p_my_reactions )
  {
    std::cerr << "No reaction data found." << std::endl;
    return EXIT_FAILURE;
  }

  //============================================================================
  // Get a view with neutrino cross sections.
  //============================================================================

  s_xpath =
    "[reactant[contains(.,'neutrino')] and .//@name = \'" +
    std::string( S_LOG10_XSEC ) +
    "\']";

  p_view =
    Libnucnet__ReacView__new(
      p_my_reactions,
      s_xpath.c_str()
    );

  //============================================================================
  // Print cross sections.
  //============================================================================

  double d_tnu = atof( argv[2] );

  fprintf( stdout, "\n\t\t\t\t\tReaction\t\t\tCross Section (cm^2)\n" );
  fprintf(
    stdout,
    "===================================================================\t====================\n"
  );

  Libnucnet__Reac__setReactionCompareFunction(
    Libnucnet__ReacView__getReac( p_view ),
    (Libnucnet__Reaction__compare_function) nnt::compare_reactions_by_string
  );

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list( Libnucnet__ReacView__getReac( p_view ) );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    fprintf(
      stdout,
      "%-67s\t%12.6e\n",
      Libnucnet__Reaction__getString( reaction.getNucnetReaction() ),
      user::compute_neutrino_cross_section(
        reaction.getNucnetReaction(),
        d_tnu
      )
    );

  }

  std::cout << std::endl;

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__ReacView__free( p_view );
  Libnucnet__Reac__free( p_my_reactions );

  return EXIT_SUCCESS;

}
