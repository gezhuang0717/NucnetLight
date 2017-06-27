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
//! \brief Example code to print out properties for a given step in the zones.
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
    std::cout << argv[0] << " input.h5 \"Step 00001\" \"exposure, n\"" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc < 4 )
  {
    std::cout << std::endl;
    std::cout << "Purpose: " << argv[0] <<
      " prints out selected properties for zones in a given step ";
    std::cout <<
      " from an hdf5 file." << std::endl;
    fprintf(
      stderr,
      "\nUsage: %s hdf5_file step property ...\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  hdf5_file = input hdf5 file\n\n"
    );
    fprintf(
      stderr,
      "  step = step to extract zone properties from\n\n"
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

  for( int i = 3; i < argc; i++ )
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

  my_zone_labels = user::hdf5::get_zone_labels( argv[1], argv[2] );

  for( size_t i = 0; i < my_zone_labels.size(); i++ )
  {

    user::hdf5::zone_labels_bimap::index<user::to>::type::iterator it =
      my_zone_labels.get<user::to>().find( i );

    std::cout <<
      it->first.tuple.get<0>().c_str() <<
      "   " <<
      it->first.tuple.get<1>().c_str() <<
      "   " <<
      it->first.tuple.get<2>().c_str();

    my_properties_hash =
      user::hdf5::get_zone_properties(
        argv[1],
        argv[2],
        it->first.tuple.get<0>().c_str(),
        it->first.tuple.get<1>().c_str(),
        it->first.tuple.get<2>().c_str()
      );

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
