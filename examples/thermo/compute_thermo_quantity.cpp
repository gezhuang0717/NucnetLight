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
//! \brief Example code to compute thermodynamic quantities from a network xml
//!        file.
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <sstream>
#include <vector>
#include "nnt/auxiliary.h"
#include "nnt/iter.h"

#include "user/thermo.h"

//##############################################################################
// check_input().
//##############################################################################

void
check_input( int argc, char **argv )
{

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {
    std::cout << std::endl;
    std::cout << argv[0] << " my_output.xml \"baryon entropy per nucleon\" \"electron internal energy density\" \"photon pressure\"" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc <= 2 ) {
    std::cout << std::endl;
    std::cout << "Purpose: " << argv[0] <<
      " computes thermodynamic quantities from a network output xml file" <<
      std::endl;

    fprintf(
      stderr, "\nUsage: %s file quantity ...\n\n", argv[0]
    );
    fprintf(
      stderr, "  file = input xml filename\n\n"
    );
    fprintf(
      stderr, "  quantity = quantity to compute (as many as desired)\n\n"
    );
    std::cout << "For an example usage, type " << std::endl << std::endl;
    std::cout << argv[0] << " --example" << std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

}

//##############################################################################
// properties routines.
//##############################################################################

std::vector<std::string> &split(
  const std::string &s,
  char delim,
  std::vector<std::string> &elems
)
{
    std::stringstream ss(s);
    std::string item;
    std::getline( ss, item, delim );
    elems.push_back( item );
    std::getline(ss, item, '\n' );
    elems.push_back( item );
    return elems;
}

std::vector<std::string> split(
  const std::string &s,
  char delim
)
{
    std::vector<std::string> elems;
    return split(s, delim, elems);
}

//##############################################################################
// main().
//##############################################################################

int
main( int argc, char * argv[] ) {

  std::vector<nnt::Zone> zones;
  std::list<std::string> quantity_list;
  Libnucnet * p_my_nucnet;

  //============================================================================
  // Check input.
  //============================================================================
  
  check_input( argc, argv ); 

  //============================================================================
  // Read input data.
  //============================================================================

  p_my_nucnet =
    Libnucnet__new_from_xml(
      argv[1],
      NULL,
      NULL,
      NULL
    );

  //============================================================================
  // Set quantities.
  //============================================================================

  for( int i = 2; i < argc; i++ )
    quantity_list.push_back( argv[i] );

  //============================================================================
  // Iterate zones.
  //============================================================================

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  BOOST_FOREACH( nnt::Zone zone, nnt::make_zone_list( p_my_nucnet ) )
  {
    zones.push_back( zone );
  }

  std::vector<std::vector<double> > props( zones.size() );
  
  for( size_t i = 0; i < zones.size(); i++ )
  {
    BOOST_FOREACH( std::string s_quantity, quantity_list )
    {
      std::vector<std::string> s_x = split( s_quantity, ' ' );

      props[i].push_back(
        user::compute_thermo_quantity( zones[i], s_x[1], s_x[0] )
      );
    }
  }

  for( size_t i = 0; i < zones.size(); i++ )
  {
    std::cout << Libnucnet__Zone__getLabel( zones[i].getNucnetZone(), 1 );

    BOOST_FOREACH( double val, props[i] )
    {
      std::cout << "  " << val;
    }

    std::cout << std::endl;

  }

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
