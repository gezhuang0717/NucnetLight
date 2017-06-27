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
//! \brief Example code to print mass fractions of all species for a given
//!        zone in a given step.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <string>
#include <iostream>

#include <boost/format.hpp>

#include "user/hdf5_routines.h"

//##############################################################################
// check_input().
//##############################################################################

void
check_input( int argc, char **argv )
{

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {
    std::cout << std::endl;
    std::cout << argv[0] << " input.h5 \"Step 00025\" 330 0 0" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc != 6 )
  {
    std::cout << std::endl;
    std::cout << "Purpose: " << argv[0] <<
      " print mass fractions for zone in a given step." << std::endl;
    fprintf(
      stderr,
      "\nUsage: %s hdf5_file step label1 label2 label3\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  hdf5_file = input hdf5 file\n\n"
    );
    fprintf(
      stderr,
      "  step = step label\n\n"
    );
    fprintf(
      stderr,
      "  label1 = first zone label\n\n"
    );
    fprintf(
      stderr,
      "  label2 = second zone label\n\n"
    );
    fprintf(
      stderr,
      "  label3 = third zone label\n\n"
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

  user::hdf5::nuclide_map my_nuclides;
  user::hdf5::zone_labels_bimap my_zone_labels;
  
  check_input( argc, argv );

  my_nuclides = user::hdf5::get_network( argv[1] );

  std::vector<user::hdf5::nuclide> my_nuclide_vector( my_nuclides.size() );

  BOOST_FOREACH( user::hdf5::nuclide_map_entry entry, my_nuclides )
  {
    my_nuclide_vector[(size_t) entry.second.iIndex] = entry.second;
  }

  user::hdf5::mass_fraction_array_t x =
    user::hdf5::get_step_mass_fractions(
      argv[1],
      argv[2]
    );

  my_zone_labels = user::hdf5::get_zone_labels( argv[1], argv[2] );

  user::hdf5::zone_labels_bimap::index<user::from>::type::iterator it =
    my_zone_labels.get<user::from>().find(
      user::hdf5::zone_labels_tuple(
        boost::make_tuple( argv[3], argv[4], argv[5] )
      )
    );

  if( it == my_zone_labels.get<user::from>().end() )
  {
    std::cerr << "Zone not found." << std::endl;
    exit( EXIT_FAILURE );
  }

  BOOST_FOREACH( user::hdf5::nuclide nuclide, my_nuclide_vector )
  {
    double d_x = x[it->second][(size_t) nuclide.iIndex];
    if( d_x > 0. )
    {
      std::cout <<
        boost::format( "%5s %3d %3d %3d   %.8e\n" ) %
        nuclide.sName %
        nuclide.iIndex %
        nuclide.iZ %
        nuclide.iA %
        d_x;
    }
  }

  return EXIT_SUCCESS;

}
