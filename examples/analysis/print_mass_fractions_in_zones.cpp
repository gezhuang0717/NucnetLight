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
//! \brief Example code to print out mass fractions in zones. 
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
// check_input().
//##############################################################################

void
check_input( int argc, char **argv )
{

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {
    std::cout << std::endl;
    std::cout << argv[0] << " my_output.xml n h1 he4" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc < 3 )
  {
    std::cout << std::endl;
    std::cout << "Purpose: " << argv[0] <<
      " prints out selected mass fractions from" << std::endl;
    std::cout <<
      "a network xml file for all zones." << std::endl;
    fprintf(
      stderr,
      "\nUsage: %s xml_file species ...\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  xml_file = network xml output file\n\n"
    );
    fprintf(
      stderr,
      "  species = name of species (enter as many as desired)\n\n"
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
  Libnucnet__Species * p_species;
  int i;
  std::list<std::string> species_list;
  
  check_input( argc, argv );

  for( i = 2; i < argc; i++ )
    species_list.push_back( argv[i] );

  p_my_nucnet = Libnucnet__new();

  Libnucnet__Nuc__updateFromXml(
    Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
    argv[1],
    NULL
  );

  Libnucnet__assignZoneDataFromXml(
    p_my_nucnet,
    argv[1],
    NULL
  );

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function)
      nnt::zone_compare_by_first_label
  );

  BOOST_FOREACH( nnt::Zone zone, nnt::make_zone_list( p_my_nucnet ) )
  {

    std::cout << Libnucnet__Zone__getLabel( zone.getNucnetZone(), 1 );

    BOOST_FOREACH( std::string s_species, species_list )
    {

      p_species =
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
          s_species.c_str()
        );

      if( !p_species )
      {
        std::cerr << std::endl;
        std::cerr << s_species << " is not a valid species!" << std::endl;
        return EXIT_FAILURE;
      }

      std::cout << "  " <<
         Libnucnet__Species__getA( p_species ) *
         Libnucnet__Zone__getSpeciesAbundance(
           zone.getNucnetZone(),
           p_species
         );

    }

    std::cout << std::endl;

  }

  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
