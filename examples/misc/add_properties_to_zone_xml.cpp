////////////////////////////////////////////////////////////////////////////////
// This file was originally written by Tianhong Yu and Bradley S. Meyer.
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
//! \brief Example code to add properties to a zone xml file.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <fstream>
#include <iostream>

#include <boost/algorithm/string.hpp>

#include <Libnucnet.h>

//##############################################################################
// check_input().
//##############################################################################

void
check_input( int argc, char **argv )
{

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {
    std::cout << std::endl;
    std::cout << argv[0] <<
      " ../../data_pub/my_net.xml ../../data_pub/zone.xml 0 0 0 prop.txt" <<
      " zone.xml" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc != 8 )
  {
    fprintf(
      stderr,
      "\nUsage: %s net_zml zone_in label_1 label_2 label_3 zone_ascii zone_out\n\n",
      argv[0]
    );
    fprintf(
      stderr, "  net_xml = input network xml filename\n\n"
    );
    fprintf(
      stderr, "  zone_in = input zone xml filename\n\n"
    );
    fprintf(
      stderr, "  label_1 = first label of desired zone\n\n"
    );
    fprintf(
      stderr, "  label_2 = second label of desired zone\n\n"
    );
    fprintf(
      stderr, "  label_3 = third label of desired zone\n\n"
    );
    fprintf(
      stderr, "  prop_ascii = input property ascii filename\n\n"
    );
    fprintf(
      stderr, "  zone_out = output zone xml filename\n\n"
    );
    std::cout << "For an example usage, type " << std::endl << std::endl;
    std::cout << argv[0] << " --example" << std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

}

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  std::ifstream prop_file;
  Libnucnet * p_my_nucnet;
  Libnucnet__Zone * p_zone;
  std::string str, s_name, s_tag1, s_tag2, s_value;
  int i;

  //============================================================================
  // Check input.
  //============================================================================

  check_input( argc, argv );

  //============================================================================
  // Read and store input.
  //============================================================================

  p_my_nucnet = Libnucnet__new();

  Libnucnet__Nuc__updateFromXml(
    Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
    argv[1],
    NULL
  );

  Libnucnet__assignZoneDataFromXml( p_my_nucnet, argv[2], NULL );

  //============================================================================
  // Get the zone.
  //============================================================================

  p_zone =
    Libnucnet__getZoneByLabels(
      p_my_nucnet,
      argv[3],
      argv[4],
      argv[5]
    );

  if( !p_zone )
  {
    std::cerr << "Zone with these labels not found.\n" << std::endl;
    return EXIT_FAILURE;
  }

  //============================================================================
  // Open property file.
  //============================================================================

  prop_file.open( argv[6] );

  if( !prop_file.is_open() )
  {
    std::cerr << "Couldn't open file " << argv[6] << "." << std::endl;
    return EXIT_FAILURE;
  }

  while( !prop_file.eof() )
  {

    s_name = "\0"; s_tag1 = "\0"; s_tag2 = "\0"; s_value = "\0";

    for( i = 0; i < 5; i++ )
    {

      if( !std::getline( prop_file, str ) ) break;

      boost::trim( str );

      if( str.empty() ) break;

      if( i == 0 ) s_name = str;
  
      if( i == 1 ) s_value = str;

      if( i == 2 )
      {
	s_tag1 = s_value;
	s_value = str;
      }

      if( i == 3 )
      {
	s_tag2 = s_value;
	s_value = str;
      }

    }

    if( !s_name.empty() )
    {

      if( s_tag1.empty() )
	Libnucnet__Zone__updateProperty(
	  p_zone,
	  s_name.c_str(),
	  NULL,
	  NULL,
	  s_value.c_str()
	);
      else
      {
	if( s_tag2.empty() )
	  Libnucnet__Zone__updateProperty(
	    p_zone,
	    s_name.c_str(),
	    s_tag1.c_str(),
	    NULL,
	    s_value.c_str()
	  );
	else
	  Libnucnet__Zone__updateProperty(
	    p_zone,
	    s_name.c_str(),
	    s_tag1.c_str(),
	    s_tag2.c_str(),
	    s_value.c_str()
	  );
      }

    }

  }

  //============================================================================
  // Write to output. 
  //============================================================================

  Libnucnet__updateZoneXmlMassFractionFormat(
    p_my_nucnet,
    "%.15e"
  );
        
  Libnucnet__writeZoneDataToXmlFile(
    p_my_nucnet,
    argv[7]
  );

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__free( p_my_nucnet );
  return EXIT_SUCCESS;

}
