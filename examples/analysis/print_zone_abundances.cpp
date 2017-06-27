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
//! \brief Example code to print out all abundances and mass fractions in
//!    zones selected by XPath.  Output is Z, A, abundance, mass fraction.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <string>
#include <iostream>
#include <Libnucnet.h>

#include "nnt/auxiliary.h"
#include "nnt/iter.h"

//##############################################################################
// Includes.
//##############################################################################

#define D_MIN    1.e-20        /* Minimum abundance to print out */

//##############################################################################
// check_input().
//##############################################################################

void
check_input( int argc, char **argv )
{

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {
    std::cout << std::endl;
    std::cout << argv[0] << " my_output.xml \"[position() >= last() - 5]\"" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc != 3 )
  {

    std::cout << std::endl;
    std::cout << "Purpose: " << argv[0] <<
      " prints out abundances from a network xml" << std::endl;
    std::cout <<
      " file for selected zones." << std::endl;
    fprintf(
      stderr,
      "\nUsage: %s xml_file zone_xpath nuc_xpath\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  xml_file = network xml output file\n\n"
    );
    fprintf(
      stderr,
      "  zone_xpath = XPath expression to select zone\n\n"
    );
    std::cout << "For an example usage, type " << std::endl << std::endl;
    std::cout << argv[0] << " --example" << std::endl << std::endl;
    exit( EXIT_FAILURE );

  }

}

//##############################################################################
// main().
//##############################################################################

int
main( int argc, char **argv )
{

  Libnucnet *p_my_nucnet;
  double d_abund;  

  check_input( argc, argv );

  p_my_nucnet = Libnucnet__new();

  Libnucnet__Nuc__updateFromXml(
    Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
    argv[1],
    NULL
  );

  Libnucnet__assignZoneDataFromXml(
    p_my_nucnet,
    argv[1],
    argv[2]
  );

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function)
      nnt::zone_compare_by_first_label
  );

  nnt::species_list_t species_list =
    nnt::make_species_list(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) )
    );

  BOOST_FOREACH( nnt::Zone zone, nnt::make_zone_list( p_my_nucnet ) )
  {

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

    BOOST_FOREACH( nnt::Species species, species_list )
    {

      d_abund = 
         Libnucnet__Zone__getSpeciesAbundance(
           zone.getNucnetZone(),
           species.getNucnetSpecies()
         );

      if( d_abund > D_MIN )
      {
        fprintf(
          stdout,
          "%5d  %5d  %.4e  %.4e\n",
          Libnucnet__Species__getZ( species.getNucnetSpecies() ),
          Libnucnet__Species__getA( species.getNucnetSpecies() ),
          d_abund,
          Libnucnet__Species__getA( species.getNucnetSpecies() ) * d_abund
        );
      }

    }

    std::cout << std::endl <<
      "1 - Xsum = " <<
      1. - Libnucnet__Zone__computeAMoment( zone.getNucnetZone(), 1 ) <<
      std::endl <<
      "Ye = " <<
      Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 ) <<
      std::endl;

    std::cout << std::endl;

  }

  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
