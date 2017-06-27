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
//! \brief Example code to create a multi-zone zone xml file from a text file.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <fstream>

#include <boost/assign/list_of.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#include <Libnucnet.h>
#include "nnt/auxiliary.h"

//##############################################################################
// check_input().
//##############################################################################

void
check_input( int argc, char **argv )
{

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {
    std::cout << std::endl;
    std::cout << argv[0]
      <<
      " ../../data_pub/my_net.xml multi_zone.txt ../../data_pub/multi_zone.xml"
      <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc != 4 )
  {
    fprintf(
      stderr,
      "\nUsage: %s nuc_file zone_ascii zone_out \n\n",
      argv[0]
    );
    fprintf(
      stderr, "  nuc_file = nuclear data xml filename\n\n"
    );
    fprintf(
      stderr, "  zone_ascii = input single zone ascii filename\n\n"
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

  typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

  Libnucnet * p_my_nucnet;
  Libnucnet__Zone * p_zone;
  std::ifstream my_file;
  std::string line, str, s_name, s_tag1, s_tag2, s_value;
  size_t i_props;

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

  //============================================================================
  // Open file.
  //============================================================================

  my_file.open( argv[2] );

  if( !my_file.is_open() )
  {
    std::cerr << "Couldn't open file." << std::endl;
    return EXIT_FAILURE;
  }

  //============================================================================
  // Read zones.
  //============================================================================

  while( std::getline( my_file, line ) )
  {

    std::vector<std::string> zone_labels =
      boost::assign::list_of("0")("0")("0");

    boost::char_separator<char> sep(",");
    tokenizer tok( line, sep);
    int j = 0;
    for( tokenizer::iterator it = tok.begin(); it != tok.end(); ++it )
    {
      std::string s_tmp = *it;
      boost::algorithm::trim( s_tmp );
      zone_labels[j] = s_tmp;
      j++;
    }

    if( j > 3 )
    {
      std::cerr << "Improper zone labels." << std::endl;
      return EXIT_FAILURE;
    }

    p_zone =
      Libnucnet__Zone__new(
        Libnucnet__getNet( p_my_nucnet ),
        zone_labels[0].c_str(),
        zone_labels[1].c_str(),
        zone_labels[2].c_str()
      );

    Libnucnet__addZone( p_my_nucnet, p_zone );

    std::getline( my_file, str );

    i_props = boost::lexical_cast<size_t>( str );

    for( size_t i = 0; i < i_props; i++ )
    {

      s_name = "\0"; s_tag1 = "\0"; s_tag2 = "\0"; s_value = "\0";

      for( size_t j = 0; j < 5; j++ )
      {

        std::getline( my_file, str );

	if( str.empty() ) break;

	if( j == 0 ) s_name = str;
    
	if( j == 1 ) s_value = str;

	if( j == 2 )
	{
	  s_tag1 = s_value;
	  s_value = str;
	}

	if( j == 3 )
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

    while( std::getline( my_file, str ) )
    {

      if( str.empty() ) break;

      boost::char_separator<char> sep1(" ");
      tokenizer tok1( str, sep1 );

      tokenizer::iterator it1 = tok1.begin();

      std::string s_tmp1 = *(it1++);
      boost::algorithm::trim( s_tmp1 );

      std::string s_tmp2 = *it1;
      boost::algorithm::trim( s_tmp2 );

      Libnucnet__Species * p_species =
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
          s_tmp1.c_str()
        );

      Libnucnet__Zone__updateSpeciesAbundance(
        p_zone,
        p_species,
        boost::lexical_cast<double>( s_tmp2 ) /
          (double) Libnucnet__Species__getA( p_species )
      );
       
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
    argv[3]
  );

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__free( p_my_nucnet );
  return EXIT_SUCCESS;

}
