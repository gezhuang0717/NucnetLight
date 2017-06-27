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
//! \brief Example code to compute Yedot in zones in a output network xml file.
//!
////////////////////////////////////////////////////////////////////////////////

#include <Libnucnet.h>

#include "nnt/auxiliary.h"
#include "nnt/string_defs.h"
#include "nnt/iter.h"

#include "user/user_rate_functions.h"
#include "user/weak_utilities.h"
#include "user/flow_utilities.h"

int main( int argc, char * argv[] ) {

  Libnucnet *p_my_nucnet;
  double d_ye;

  boost::tuple<double, double, double, double, double> t;

  //============================================================================
  // Check input.
  //============================================================================

   if ( argc != 3 ) {
      fprintf(
        stderr, "\nUsage: %s xml_filename zone_xpath \n\n", argv[0]
      );
      fprintf(
        stderr, "  xml_filename = input xml filename\n\n"
      );
      fprintf(
        stderr, "  zone_xpath = XPATH to select zones\n\n"
      );
      return EXIT_FAILURE;
   }

  //============================================================================
  // Read file and exit if not present.
  //============================================================================

  p_my_nucnet = Libnucnet__new_from_xml( argv[1], NULL, NULL, argv[2] );

  if( !p_my_nucnet ) {
    fprintf( stderr, "Input data not read!\n" );
    return EXIT_FAILURE;
  }

  //============================================================================
  // Register user rate functions.
  //============================================================================

  user::register_my_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

  //============================================================================
  // Set the zone views.  This greatly speeds up the calculation.
  //============================================================================

  user::set_weak_views_in_zones( p_my_nucnet );

  //============================================================================
  // Iterate the zones.
  //============================================================================
  
  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  nnt::zone_list_t zone_list = nnt::make_zone_list( p_my_nucnet );

  BOOST_FOREACH( nnt::Zone zone, zone_list )
  {

    if( !zone.hasProperty( nnt::s_MU_NUE_KT ) )
    {
      zone.updateProperty(
        nnt::s_MU_NUE_KT,
        boost::lexical_cast<std::string>( GSL_NEGINF )
      );
    }

    d_ye = 
      Libnucnet__Zone__computeZMoment(
	zone.getNucnetZone(), 1
      );

    user::update_my_rate_functions_data( zone );

    t = user::compute_all_yedot( d_ye, &zone );

    fprintf(
      stdout,
      "\n ye = %g\n yedot =  %g\n tau_ye =  %g\n\n",
      d_ye,
      t.get<4>(),
      d_ye / t.get<4>()
    );

    fprintf(
      stdout,
      "\
       Beta-minus yedot       = %g\n\
       Beta-plus yedot        = %g\n\
       Electron-capture yedot = %g\n\
       Positron-capture yedot = %g\n\
      ",
      t.get<0>(),
      t.get<1>(),
      t.get<2>(),
      t.get<3>()
    );

    std::cout << std::endl;

  }

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
