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
//! \brief Example code to print out properties in a network xml file.
////////////////////////////////////////////////////////////////////////////////

/*##############################################################################
// Includes.
//############################################################################*/

#include <iostream>
#include <string>

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>

#include "nnt/auxiliary.h"
#include "nnt/iter.h"

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
    std::cout << argv[0] << " my_output.xml time t9 rho \"exposure, n\"" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc < 3 )
  {
    std::cout << std::endl;
    std::cout << "Purpose: " << argv[0] <<
      " prints out selected properties from a network xml" << std::endl;
    std::cout <<
      " file for all zones." << std::endl;
    fprintf(
      stderr,
      "\nUsage: %s xml_file prop ...\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  xml_file = network xml output file\n\n"
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
// count_properties().
//##############################################################################

void
count_properties(
  const char * s_name,
  const char * s_tag1,
  const char * s_tag2,
  const char * s_value,
  int &i
)
{

  i++;

} 

//##############################################################################
// count_instances().
//##############################################################################

int
count_instances( nnt::Zone& zone, std::string property1, std::string property2 )
{

  int i = 0;

  Libnucnet__Zone__iterateOptionalProperties(
    zone.getNucnetZone(),
    property1.c_str(),
    property2.c_str(),
    NULL,
    (Libnucnet__Zone__optional_property_iterate_function) count_properties,
    &i
  );

  return i;

}

int
count_instances( nnt::Zone& zone, std::string property )
{

  int i = 0;

  Libnucnet__Zone__iterateOptionalProperties(
    zone.getNucnetZone(),
    property.c_str(),
    NULL,
    NULL,
    (Libnucnet__Zone__optional_property_iterate_function) count_properties,
    &i
  );

  return i;

}

//##############################################################################
// main().
//##############################################################################

int
main( int argc, char **argv )
{

  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

  Libnucnet *p_my_nucnet;
  std::list<std::vector<std::string> > property_list;

  check_input( argc, argv );

  for( int i = 2; i < argc; i++ )
  {

    std::string s_arg = argv[i];
    std::vector<std::string> prop_tags;

    boost::char_separator<char> sep(",");
    tokenizer tok( s_arg, sep);
    int j = 0;
    for (tokenizer::iterator it = tok.begin(); it != tok.end(); ++it)
    {
      std::string s_tmp = *it;
      boost::algorithm::trim( s_tmp );
      prop_tags.push_back( s_tmp );
      j++;
    }

    if( j > 3 )
    {
      std::cerr << "Invalid property and tags." << std::endl;
      return EXIT_FAILURE;
    }

    property_list.push_back( prop_tags );

  }

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
      zone_compare_by_first_label
  );

  BOOST_FOREACH( nnt::Zone zone, nnt::make_zone_list( p_my_nucnet ) )
  {

    std::cout << Libnucnet__Zone__getLabel( zone.getNucnetZone(), 1 );

    BOOST_FOREACH( std::vector<std::string> property, property_list )
    {

      if( property.size() == 1 )
      {
        if( count_instances( zone, property[0] ) > 1 )
        {
          std::cerr << "More than one instance of property " <<
            property[0] << std::endl;
          return EXIT_FAILURE;
        }
        std::cout << "  " << zone.getProperty( property[0] );
      }

      if( property.size() == 2 )
      {
        if( count_instances( zone, property[0], property[1] ) > 1 )
        {
          std::cerr << "More than one instance of property " <<
            property[0] << "  " << property[1] << std::endl;
          return EXIT_FAILURE;
        }
        std::cout << "  " <<
          zone.getProperty( property[0], property[1] );
      }

      if( property.size() == 3 )
      {
        std::cout << "  " <<
          zone.getProperty( property[0], property[1], property[2] );
      }

    }

    std::cout << std::endl;

  }

  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
