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
//! \brief Example code to print out abundance moments in zones. 
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <string>
#include <iostream>
#include <Libnucnet.h>

#include <boost/format.hpp>

#include "nnt/auxiliary.h"
#include "nnt/iter.h"
#include "user/network_utilities.h"

//##############################################################################
// check_input().
//##############################################################################

void
check_input( int argc, char **argv )
{

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {
    std::cout << std::endl;
    std::cout << argv[0] << " my_output.xml \"\" z 1" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc < 5 || argc > 6 )
  {
    std::cout << std::endl;
    std::cout << "Purpose: " << argv[0] <<
      " prints out selected mass fractions from" << std::endl;
    std::cout <<
      "a network xml file for all zones." << std::endl;
    fprintf(
      stderr,
      "\nUsage: %s xml_file nuc_xpath nucleon exponent zone_xpath\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  xml_file = network xml output file\n\n"
    );
    fprintf(
      stderr,
      "  nuc_xpath = XPath expression to select nuclei\n\n"
    );
    fprintf(
      stderr,
      "  nucleon = z, n, or a for moment\n\n"
    );
    fprintf(
      stderr,
      "  exponent = exponent on z, n, or a for moment\n\n"
    );
    fprintf(
      stderr,
      "  zone_xpath = XPath expression to select zones (optional)\n\n"
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
  
  check_input( argc, argv );

  p_my_nucnet = Libnucnet__new();

  Libnucnet__Nuc__updateFromXml(
    Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
    argv[1],
    NULL
  );

  if( argc == 5 )
  {
    Libnucnet__assignZoneDataFromXml(
      p_my_nucnet,
      argv[1],
      NULL
    );
  }
  else
  {
    Libnucnet__assignZoneDataFromXml(
      p_my_nucnet,
      argv[1],
      argv[5]
    );
  }

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function)
      nnt::zone_compare_by_first_label
  );

  nnt::zone_list_t zone_list = nnt::make_zone_list( p_my_nucnet );

  boost::format fmt( "%d  %g\n" );

  BOOST_FOREACH( nnt::Zone zone, zone_list )
  {

    if( !( zone == *(zone_list.begin()) ) )
    {
      Libnucnet__Zone__copy_net_views(
        zone.getNucnetZone(),
        (*zone_list.begin()).getNucnetZone()
      );
    }

    fmt %
      Libnucnet__Zone__getLabel( zone.getNucnetZone(), 1 ) %
      user::compute_cluster_abundance_moment(
        zone,
        argv[2],
        argv[3],
        atoi( argv[4] )
      );

    std::cout << fmt.str();
        
  }

  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
