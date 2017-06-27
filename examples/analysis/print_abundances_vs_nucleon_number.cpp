////////////////////////////////////////////////////////////////////////////////
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
//! \brief Example code to print out abundances vs. nucleon number
//!        in a network xml file.
////////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <Libnucnet.h>

#include "nnt/auxiliary.h"
#include "nnt/iter.h"

int
main( int argc, char * argv[] ) {

  Libnucnet * p_my_nucnet;
  gsl_vector * p_vector;
  size_t i;

  //============================================================================
  // Check input.
  //============================================================================

  if ( argc!= 4 ) {
      fprintf(
        stderr, "\nUsage: %s file nucleon zone_xpath\n\n", argv[0]
      );
      fprintf(
        stderr, "  file = input xml filename\n\n"
      );
      fprintf(
        stderr, "  nucleon = nucleon type (z, n, or a)\n\n"
      );
      fprintf(
        stderr, "  zone_xpath = XPath to select zone\n\n"
      );

      return EXIT_FAILURE;
  }

  //============================================================================
  // Read input data.
  //============================================================================

  p_my_nucnet = Libnucnet__new_from_xml( argv[1], NULL, NULL, argv[3] );

  //============================================================================
  // Print out largest nucleon number.
  //============================================================================

  std::cout <<
    Libnucnet__Nuc__getLargestNucleonNumber(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
      argv[2]
    ) <<
    std::endl;

  //============================================================================
  // Iterate zones.
  //============================================================================

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  BOOST_FOREACH( nnt::Zone zone, nnt::make_zone_list( p_my_nucnet ) )
  {

    p_vector =
      Libnucnet__Zone__getSummedAbundances(
        zone.getNucnetZone(),
        argv[2]
      );

    if( zone.hasProperty( nnt::s_TIME ) )
      std::cout <<
        "time(s) = " << zone.getProperty( nnt::s_TIME ) << " " <<
        "t9 = " << zone.getProperty( nnt::s_T9 ) << " " <<
        "rho(g/cc) = " << zone.getProperty( nnt::s_RHO ) << " " <<
        std::endl;
    else
      std::cout <<
        "t9 = " << zone.getProperty( nnt::s_T9 ) << " " <<
        "rho(g/cc) = " << zone.getProperty( nnt::s_RHO ) << " " <<
        std::endl;

    for( i = 0; i < WnMatrix__get_gsl_vector_size( p_vector ); i++ )
    {
      fprintf(
	stdout,
	"%lu  %.4e\n",
	(unsigned long) i,
	gsl_vector_get( p_vector, i )
      );
    }

    gsl_vector_free( p_vector );

    std::cout << std::endl;

  }

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}

