//////////////////////////////////////////////////////////////////////////////
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
//! \brief Example code to compute the rate of change of the entropy per
//!    nucleon in zones in a output network xml file.
//!
////////////////////////////////////////////////////////////////////////////////

#include <omp.h>

#include <Libnucnet.h>

#include "nnt/auxiliary.h"
#include "nnt/string_defs.h"
#include "nnt/iter.h"

#include "user/aa522a25.h"
#include "user/flow_utilities.h"
#include "user/user_rate_functions.h"
#include "user/thermo.h"

#define APPROXIMATE_WEAK_RATES  1

//##############################################################################
// check_input().
//##############################################################################

void
check_input( int argc, char **argv )
{

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {
    std::cout << std::endl;
    std::cout <<
      argv[0] <<
      " my_output.xml \"[optional_properties/property[@name = 't9'] > 1]\"" <<
    std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc != 3 )
  {
    std::cout << std::endl;
    std::cout << "Purpose: " << argv[0] <<
      " computes entropy generation rate per nucleon from a network output xml file" <<
      std::endl;
    fprintf(
      stderr, "\nUsage: %s xml_filename zone_xpath \n\n", argv[0]
    );
    fprintf(
      stderr, "  xml_filename = input xml filename\n\n"
    );
    fprintf(
      stderr, "  zone_xpath = XPATH to select zones\n\n"
    );
    exit( EXIT_FAILURE );
  }

}

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  Libnucnet *p_my_nucnet;
  Libnucnet__NetView * p_view;
  std::vector<nnt::Zone> zones;

  //============================================================================
  // Check input.
  //============================================================================

  check_input( argc, argv );

  //============================================================================
  // Read file and exit if not present.
  //============================================================================

  p_my_nucnet = Libnucnet__new_from_xml( argv[1], NULL, NULL, argv[2] );

  if( !p_my_nucnet ) {
    fprintf( stderr, "Input data not read!\n" );
    return EXIT_FAILURE;
  }

  //============================================================================
  // Update with approximate rates.
  //============================================================================

  if( APPROXIMATE_WEAK_RATES )
    user::aa522a25__update_net(
      Libnucnet__getNet( p_my_nucnet )
    );

  //============================================================================
  // Register user rate functions.
  //============================================================================

  user::register_my_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

  //============================================================================
  // Get the view.
  //============================================================================
  
  p_view =
    Libnucnet__NetView__new(
      Libnucnet__getNet( p_my_nucnet ),
      "",
      ""
    );

  //============================================================================
  // Iterate the zones.
  //============================================================================
  
  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  BOOST_FOREACH( nnt::Zone zone, nnt::make_zone_list( p_my_nucnet ) )
  {
    zones.push_back( zone );
  }

  #pragma omp parallel for schedule( dynamic, 1 )
  for( size_t i = 0; i < zones.size(); i++ )
  {

    Libnucnet__Zone__updateNetView(
      zones[i].getNucnetZone(),
      "",
      "",
      NULL,
      Libnucnet__NetView__copy( p_view )
    );

    if(
      !Libnucnet__Zone__getProperty(
        zones[i].getNucnetZone(),
        nnt::s_MU_NUE_KT,
        NULL,
        NULL
      )
    )
      zones[i].updateProperty(
	nnt::s_MU_NUE_KT,
	boost::lexical_cast<std::string>( GSL_NEGINF )
      );

    user::update_my_rate_functions_data( zones[i] );

    zones[i].updateProperty(
      nnt::s_SDOT,
      boost::lexical_cast<std::string>(
        user::compute_entropy_change_rate(
          zones[i],
          zones[i].getNetView( "", "" )
        )
      )
    );

  }

  BOOST_FOREACH( nnt::Zone zone, zones )
  {
    std::cout <<
      Libnucnet__Zone__getLabel( zone.getNucnetZone(), 1 ) << "  " <<
      zone.getProperty( nnt::s_SDOT ) << std::endl;
  }

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__NetView__free( p_view );
  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
