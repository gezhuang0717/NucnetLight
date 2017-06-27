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
//! \brief Example code to print out properties for a given zone in all steps.
////////////////////////////////////////////////////////////////////////////////

/*##############################################################################
// Includes.
//############################################################################*/

#include <iostream>
#include <string>

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

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
    std::cout << argv[0] << " input.h5 535 0 0 \"exposure, n\"" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc < 6 )
  {
    std::cout << std::endl;
    std::cout << "Purpose: " << argv[0] <<
      " prints out selected properties for a zone in all steps ";
    std::cout <<
      "from an hdf5 file." << std::endl;
    fprintf(
      stderr,
      "\nUsage: %s hdf5_file label1 label2 label3 property ...\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  hdf5_file = input hdf5 file\n\n"
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
    fprintf(
      stderr,
      "  prop = name of property (enter as many as desired)\n\n"
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

  user::hdf5::zone_labels_bimap my_zone_labels;
  user::hdf5::zone_properties_hash my_properties_hash;
  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

  std::list<std::vector<std::string> > property_list;

  check_input( argc, argv );

  for( int i = 5; i < argc; i++ )
  {

    std::string s_arg = argv[i];
    std::vector<std::string> prop_tags;

    for( size_t j = 0; j < 3; j++ ) prop_tags.push_back( "0" );

    boost::char_separator<char> sep(",");
    tokenizer tok( s_arg, sep);
    int j = 0;
    for (tokenizer::iterator it = tok.begin(); it != tok.end(); ++it)
    {
      if( j > 3 )
      {
        std::cerr << "Invalid property and tags." << std::endl;
        return EXIT_FAILURE;
      }
      std::string s_tmp = *it;
      boost::algorithm::trim( s_tmp );
      prop_tags[j] = s_tmp;
      j++;
    }

    property_list.push_back( prop_tags );

  }

  size_t i_steps = user::hdf5::count_groups_in_file( argv[1] );

  for( size_t i = 0; i < i_steps; i++ )
  {

    std::string s_g = user::hdf5::get_step_label( i );

    my_zone_labels = user::hdf5::get_zone_labels( argv[1], s_g.c_str() );

    my_properties_hash =
      user::hdf5::get_zone_properties(
        argv[1],
        s_g.c_str(),
        argv[2],
        argv[3],
        argv[4]
      );

    std::cout << i << "   ";

    BOOST_FOREACH( std::vector<std::string> property, property_list )
    {

      user::hdf5::zone_properties_hash::iterator it2 =
        my_properties_hash.find(
          boost::make_tuple(
            property[0],
            property[1],
            property[2]
          )
        ); 

      std::cout <<
        "   " <<
        it2->sValue;

    }

    std::cout << std::endl;

  }

  return EXIT_SUCCESS;

}
