////////////////////////////////////////////////////////////////////////////////
// This file was originally written by Tianhong Yu. 
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
//! \brief Example code to extract zones from xml files and output
//!    to a new xml file.
//! The zone list file should contain each xml file name in each line.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <string>
#include <iostream>
#include <fstream>

#include <boost/tuple/tuple.hpp>

#include <Libnucnet.h>

#include "nnt/auxiliary.h"
#include "nnt/iter.h"

//##############################################################################
// main().
//##############################################################################

int
main( int argc, char **argv )
{

  Libnucnet *p_my_nucnet;
  std::ifstream my_file;
  std::string s_xml_input;
  std::vector<boost::tuple<std::string,std::string,std::string> > old_labels;
  nnt::zone_list_t zone_list;

  if( argc != 4 )
  {
    fprintf(
      stderr,
      "\nUsage: %s input_file_list xpath xml_file\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  input_file_list = a text file including all input xml files\n\n"
    );
    fprintf(
      stderr,
      "  xpath = xpath to select zones\n\n"
    );
    fprintf(
      stderr,
      "  xml_file = network xml output file\n\n"
    );
    return EXIT_FAILURE;
  }

  p_my_nucnet = Libnucnet__new();

  //============================================================================
  // Open the file.  Each line of the file should have the name of an input
  // xml file.
  //============================================================================

  my_file.open( argv[1] );

  if( !my_file.is_open() )
  {
    std::cerr << "Couldn't open file." << std::endl;
    return EXIT_FAILURE;
  }

  //============================================================================
  // Loop over the files.  Skip empty lines.
  //============================================================================

  while( std::getline( my_file, s_xml_input ) )
  {

    if( !s_xml_input.empty() )
    {

      //------------------------------------------------------------------------
      // Update Nuc data with the first zone.
      //------------------------------------------------------------------------

      if(
	Libnucnet__Nuc__getNumberOfSpecies(
	  Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) )
	) == 0
      )
      {
	Libnucnet__Net__updateFromXml(
	  Libnucnet__getNet( p_my_nucnet ),
	  s_xml_input.c_str(),
	  NULL,
	  NULL
	);
      }

      //------------------------------------------------------------------------
      // Assign zones.
      //------------------------------------------------------------------------

      size_t i = Libnucnet__getNumberOfZones( p_my_nucnet );

      Libnucnet__assignZoneDataFromXml(
	p_my_nucnet,
	s_xml_input.c_str(),
	argv[2]
      );

      //------------------------------------------------------------------------
      // Relabel zones.
      //------------------------------------------------------------------------

      Libnucnet__setZoneCompareFunction(
	p_my_nucnet,
	(Libnucnet__Zone__compare_function)
	  nnt::zone_compare_by_first_label 
      );

      zone_list = nnt::make_zone_list( p_my_nucnet );

      BOOST_FOREACH( nnt::Zone zone, zone_list )
      {

	if(
	  std::string( Libnucnet__Zone__getLabel( zone.getNucnetZone(), 2 ) ) !=
	  "added"
	)
	{
	  old_labels.push_back(
	    boost::make_tuple(
	      Libnucnet__Zone__getLabel( zone.getNucnetZone(), 1 ),
	      Libnucnet__Zone__getLabel( zone.getNucnetZone(), 2 ),
	      Libnucnet__Zone__getLabel( zone.getNucnetZone(), 3 )
	    )
	  );
	  Libnucnet__relabelZone(
	    p_my_nucnet,
	    zone.getNucnetZone(),
	    boost::lexical_cast<std::string>( i++ ).c_str(),
	    "added",
	    "0"
	  );
	}

      }

    }

  }

  my_file.close();

  //============================================================================
  // Final zone relabel.  Here we relabel the "added" with the original
  // label 1.
  //============================================================================

  zone_list = nnt::make_zone_list( p_my_nucnet );

  BOOST_FOREACH( nnt::Zone zone, zone_list )
  {

    std::string s_label =
      Libnucnet__Zone__getLabel( zone.getNucnetZone(), 1 );

    Libnucnet__relabelZone(
      p_my_nucnet,
      zone.getNucnetZone(),
      s_label.c_str(),
      old_labels[
        boost::lexical_cast<size_t>(
          Libnucnet__Zone__getLabel( zone.getNucnetZone(), 1 )
        )
      ].get<0>().c_str(),
      "0"
    );

  }

  //============================================================================
  // Write output.
  //============================================================================

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function)
      nnt::zone_compare_by_first_label 
  );

  Libnucnet__writeToXmlFile( p_my_nucnet, argv[3] );

  //============================================================================
  // Clean up. 
  //============================================================================

  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
