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
    std::cout << argv[0] << " my_file.hdf5 \"Step 00001\" h1 he4" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc < 3 )
  {
    std::cout << std::endl;
    std::cout << "Purpose: " << argv[0] <<
      " prints out selected mass fractions from ";
    std::cout <<
      "a network hdf5 file for all steps for a given zone." << std::endl;
    fprintf(
      stderr,
      "\nUsage: %s hdf5_file step species ...\n\n",
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

  for( int i = 3; i < argc; i++ )
    species_list.push_back( argv[i] );

  my_nuclides = user::hdf5::get_network( argv[1] );

  my_zone_labels = user::hdf5::get_zone_labels( argv[1], argv[2] );
       
  user::hdf5::mass_fraction_array_t x =
    user::hdf5::get_step_mass_fractions(
      argv[1],
      argv[2]
    );

  for( size_t i = 0; i < my_zone_labels.size(); i++ )
  {

    user::hdf5::zone_labels_bimap::index<user::to>::type::iterator it =
      my_zone_labels.get<user::to>().find( i );

    if( it == my_zone_labels.get<user::to>().end() )
    {
      std::cerr << "Zone not found." << std::endl;
      exit( EXIT_FAILURE );
    }

    std::cout <<
      it->first.tuple.get<0>() <<
      "   " <<
      it->first.tuple.get<1>() <<
      "   " <<
      it->first.tuple.get<2>();

    BOOST_FOREACH( std::string s_nuclide, species_list )
    {

      if( my_nuclides.find( s_nuclide ) == my_nuclides.end() )
      {
        std::cerr <<
          "  " << s_nuclide << " not present in collection." << std::endl;
        exit( EXIT_FAILURE );
      }

      std::cout <<
        "   " <<
        x[i][(size_t) my_nuclides[s_nuclide].iIndex];

    }

    std::cout << std::endl;

  }
        
  return EXIT_SUCCESS;

}
