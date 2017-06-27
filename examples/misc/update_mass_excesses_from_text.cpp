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
//! \brief Example code to update mass excesses in a nuclide or network xml
//!        data file.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <fstream>
#include <iostream>

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

#include <nnt/auxiliary.h>
#include <Libnucnet.h>

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  Libnucnet__Net * p_net;
  Libnucnet__Nuc * p_nuc;
  std::string s_type, line, s_state;
  double d_mass_excess;
  typedef boost::tokenizer<boost::char_separator<char> > my_tok;
  boost::char_separator<char> sep( " " );

  //============================================================================
  // Check input.
  //============================================================================

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {
    std::cout << std::endl;
    std::cout << argv[0] <<
      " ../../data_pub/my_net.xml net nuc.txt foo new_net.xml" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc != 6 )
  {
    fprintf(
      stderr,
      "\nUsage: %s nuc_xml type nuc_text source output_xml\n\n",
      argv[0]
    );
    fprintf(
      stderr, "  nuc_xml = original input xml filename\n\n"
    );
    fprintf(
      stderr, "  type = type of file (\"nuc\" or \"net\")\n\n"
    );
    fprintf(
      stderr, "  nuc_text = text file with new data\n\n"
    );
    fprintf(
      stderr, "  source = name of the source of the mass excess data\n\n"
    );
    fprintf(
      stderr, "  output_xml = output xml file\n\n"
    );
    std::cout << "For an example usage, type " << std::endl << std::endl;
    std::cout << argv[0] << " --example" << std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  //============================================================================
  // Set type.
  //============================================================================

  s_type = argv[2];

  if( s_type != "nuc" && s_type != "net" )
  {
    std::cerr << "Invalid file type (should be \"nuc\" or \"net\")."
              << std::endl;
    return EXIT_FAILURE;
  }

  //============================================================================
  // Read and store input.
  //============================================================================

  if( s_type == "nuc" )
  {
    p_nuc = Libnucnet__Nuc__new_from_xml( argv[1], NULL );
  }
  else if( s_type == "net" )
  {
    p_net = Libnucnet__Net__new_from_xml( argv[1], NULL, NULL );
    p_nuc = Libnucnet__Net__getNuc( p_net );
  }

  //============================================================================
  // Open file.
  //============================================================================

  std::ifstream input_file( argv[3] );

  if( !input_file.good() )
  {
    std::cerr << "Couldn't open file " << argv[3] << "." << std::endl;
    return EXIT_FAILURE;
  }

  //============================================================================
  // Loop over input and update mass excesses.
  //============================================================================

  while( std::getline( input_file, line ) )
  {

    my_tok tok( line, sep );

    std::vector<std::string> items;

    for( my_tok::iterator it = tok.begin(); it != tok.end();++it )
    {
      items.push_back( *it );
    }

    if( items.size() > 2 )   // Two items or fewer is invalid input--skip line.
    {

      if( items.size() == 3 )
      {
	s_state = "";
	d_mass_excess = boost::lexical_cast<double>( items[2] );
      }
      else
      {
	s_state = items[2];
	d_mass_excess = boost::lexical_cast<double>( items[3] );
      }

      Libnucnet__Species * p_species =
	Libnucnet__Nuc__getSpeciesByZA(
	  p_nuc,
	  boost::lexical_cast<unsigned int>( items[0] ),
	  boost::lexical_cast<unsigned int>( items[1] ),
	  s_state.c_str()
	);

      if( p_species )
      {
	Libnucnet__Species__updateMassExcess( p_species, d_mass_excess );
	Libnucnet__Species__updateSource( p_species, argv[4] );
      }
      else
      {
	std::cerr << "Species: " << items[0] << ", " << items[1] << s_state <<
	  " not in original collection." << std::endl;
	return EXIT_FAILURE;
      }

    }

  }

  input_file.close();

  //============================================================================
  // Write to output.
  //============================================================================

  if( s_type == "nuc" )
    Libnucnet__Nuc__writeToXmlFile( p_nuc, argv[5] );
  else
    Libnucnet__Net__writeToXmlFile( p_net, argv[5] );

  //============================================================================
  // Clean up and exit.
  //============================================================================

  if( s_type == "nuc" )
    Libnucnet__Nuc__free( p_nuc );
  else
    Libnucnet__Net__free( p_net );

  return EXIT_SUCCESS;

}
