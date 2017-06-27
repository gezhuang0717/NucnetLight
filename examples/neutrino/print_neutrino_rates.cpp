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
//! \brief Example code to print neutrino rates.
////////////////////////////////////////////////////////////////////////////////

#include <Libnucnet__Reac.h>
#include "nnt/string_defs.h"
#include "nnt/auxiliary.h"
#include "nnt/iter.h"

#include "user/neutrino_rate_functions.h"

int main( int argc, char * argv[] ) {

  Libnucnet *p_my_nucnet;
  double d_rate;
  double d_t9 = 1.;   // Default neutrino rates do not need t9 but rather
                      // the neutrino temperature

  //============================================================================
  // Check input.
  //============================================================================

  if ( argc != 4 && argc!= 5 ) {
     fprintf(
       stderr, "\nUsage: %s neutrino_xml zone_xml r xpath\n\n", argv[0]
     );
     fprintf(
       stderr, "  neutrino_xml = input neutrino data xml filename\n\n"
     );
     fprintf(
       stderr, "  zone_xml = input zone data xml filename\n\n"
     );
     fprintf(
       stderr, "  r = input radius (cm)\n\n"
     );
     fprintf(
       stderr, "  xpath = reaction xpath expression (optional)\n\n"
     );
     return EXIT_FAILURE;
  }

  //============================================================================
  // Read file and exit if not present.
  //============================================================================

  p_my_nucnet = Libnucnet__new();

  if ( argc == 4 )
    Libnucnet__Reac__updateFromXml(
      Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) ),
      argv[1],
      NULL
    );
  else
    Libnucnet__Reac__updateFromXml(
      Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) ),
      argv[1],
      argv[4]
    );

    if( !p_my_nucnet ) exit(0);

  //============================================================================
  // Get zone data, which includes neutrino temperatures and luminosities.
  //============================================================================

  Libnucnet__assignZoneDataFromXml( p_my_nucnet, argv[2], NULL );

  if( !Libnucnet__getZoneByLabels( p_my_nucnet, "0", "0", "0" ) )
  {
    std::cerr << "No such zone." << std::endl;
    return EXIT_FAILURE;
  }

  nnt::Zone zone;

  zone.setNucnetZone(
    Libnucnet__getZoneByLabels( p_my_nucnet, "0", "0", "0" )
  );

  //============================================================================
  // Radius.
  //============================================================================

  zone.updateProperty(
    nnt::s_RADIUS,
    argv[3]
  );

  //============================================================================
  // Print rates.
  //============================================================================

  fprintf( stdout, "\n\t\t\t\t\tReaction\t\t\tRate\n" );
  fprintf(
    stdout,
    "===================================================================\t============\n"
  );

  Libnucnet__Reac__setReactionCompareFunction(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) ),
    (Libnucnet__Reaction__compare_function) nnt::compare_reactions_by_string
  );

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list( 
      Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
    );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    if(
	strcmp(
	  Libnucnet__Reaction__getRateFunctionKey( reaction.getNucnetReaction() ),
	  NU_NUCL
	) == 0
    )
    {
      d_rate =
	user::nu_nucl_function( reaction.getNucnetReaction(), d_t9, zone );
    }
    else if(
	strcmp(
	  Libnucnet__Reaction__getRateFunctionKey( reaction.getNucnetReaction() ),
	  NU_N_CAPTURE
	) == 0
    )
    {
      d_rate =
	user::nu_n_capture_function(
	  reaction.getNucnetReaction(),
	  d_t9,
	  zone
	);
    }
    else if(
	strcmp(
	  Libnucnet__Reaction__getRateFunctionKey(
            reaction.getNucnetReaction()
          ),
	  NU_P_CAPTURE
	) == 0
    )
    {
      d_rate =
	user::nu_p_capture_function(
	  reaction.getNucnetReaction(),
	  d_t9,
	  zone
	);
    }
    else
    {
      fprintf( stderr, "Rate function not found.\n" );
      exit( EXIT_FAILURE );
    }

    fprintf(
      stdout,
      "%-67s\t%12.6e\n",
      Libnucnet__Reaction__getString( reaction.getNucnetReaction() ),
      d_rate
    );

  }

  fprintf( stdout, "\n" );

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
