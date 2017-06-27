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
//! \brief Example code to print out mass fractions from a zone in the
//!        steps.
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
    std::cout << argv[0] << " my_file.hdf5 500 0 0 h1 he4" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc < 6 )
  {
    std::cout << std::endl;
    std::cout << "Purpose: " << argv[0] <<
      " prints out selected mass fractions from ";
    std::cout <<
      "a network hdf5 file for all steps for a given zone." << std::endl;
    fprintf(
      stderr,
      "\nUsage: %s hdf5_file zone_label1 zone_label2 zone_label3 species ...\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  hdf5_file = input hdf5 file\n\n"
    );
    fprintf(
      stderr,
      "  zone_label1 = first label of zone\n\n"
    );
    fprintf(
      stderr,
      "  zone_label2 = second label of zone\n\n"
    );
    fprintf(
      stderr,
      "  zone_label3 = third label of zone\n\n"
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

  user::hdf5::nuclide_map my_nuclides;
  user::hdf5::zone_labels_bimap my_zone_labels;
  std::list<std::string> species_list;
  
  check_input( argc, argv );

  for( int i = 5; i < argc; i++ )
    species_list.push_back( argv[i] );

  size_t i_steps = user::hdf5::count_groups_in_file( argv[1] );

  my_nuclides = user::hdf5::get_network( argv[1] );

  for( size_t i = 0; i < i_steps; i++ )
  {

    std::string s_step = user::hdf5::get_step_label( i );

    my_zone_labels = user::hdf5::get_zone_labels( argv[1], s_step.c_str() );
       
    user::hdf5::zone_labels_bimap::index<user::from>::type::iterator it =
      my_zone_labels.get<user::from>().find(
        user::hdf5::zone_labels_tuple(
          boost::make_tuple( argv[2], argv[3], argv[4] )
        )
      );

    if( it == my_zone_labels.get<user::from>().end() )
    {
      std::cerr << "Zone not found." << std::endl;
      exit( EXIT_FAILURE );
    }

    user::hdf5::mass_fraction_array_t x =
      user::hdf5::get_step_mass_fractions(
        argv[1],
        s_step.c_str()
      );

    BOOST_FOREACH( std::string s_nuclide, species_list )
    {

      if( my_nuclides.find( s_nuclide ) == my_nuclides.end() )
      {
        std::cerr << s_nuclide << " not present in collection." << std::endl;
        exit( EXIT_FAILURE );
      }

      std::cout <<
        boost::format( "%.4e   " ) %
        x[it->second][(size_t) my_nuclides[s_nuclide].iIndex];

    }

    std::cout << std::endl;

  }
        
  return EXIT_SUCCESS;

}
