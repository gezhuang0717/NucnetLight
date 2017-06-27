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
//! \brief Example code to print out sound speed for zones in a network xml
//!    file.
////////////////////////////////////////////////////////////////////////////////

/*##############################################################################
// Includes.
//############################################################################*/

#include <omp.h>
#include <iostream>
#include <string>

#include "nnt/auxiliary.h"
#include "nnt/iter.h"

#include "user/thermo.h"

using namespace nnt;

//##############################################################################
// check_input().
//##############################################################################

void
check_input( int argc, char **argv )
{

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {
    std::cout << std::endl;
    std::cout << argv[0] << " my_output.xml \"[position() >= last() - 10]\"" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc < 3 )
  {
    std::cout << std::endl;
    std::cout << "Purpose: " << argv[0] <<
      " computes sound speed for selected zones from a network xml file." <<
      std::endl;
    fprintf(
      stderr,
      "\nUsage: %s xml_file zone_xpath\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  xml_file = network xml output file\n\n"
    );
    fprintf(
      stderr,
      "  zone_xpath = xpath expression to select zones\n\n"
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
  std::vector<Zone> zones;

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
      zone_compare_by_first_label
  );

  BOOST_FOREACH( nnt::Zone zone, nnt::make_zone_list( p_my_nucnet ) )
  {
    zones.push_back( zone );
  }

  for( size_t i = 0; i < zones.size(); i++ )
  {
    zones[i].updateProperty(
      nnt::s_SOUND_SPEED,
      boost::lexical_cast<std::string>(
        user::compute_sound_speed( zones[i] )
      )
    );
  }

  BOOST_FOREACH( nnt::Zone zone, zones )
  {
    std::cout <<
      Libnucnet__Zone__getLabel( zone.getNucnetZone(), 1 ) << "  " <<
      zone.getProperty( nnt::s_SOUND_SPEED ) << std::endl;
  }
    
  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
