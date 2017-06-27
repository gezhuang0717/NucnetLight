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
//! \brief Example code to print out all optional properties in
//!    zones selected by XPath.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <string>
#include <iostream>
#include <Libnucnet.h>

#include "nnt/auxiliary.h"
#include "nnt/iter.h"

//##############################################################################
// Prototype.
//##############################################################################

void
print_optional_properties(
  const char *,
  const char *,
  const char *,
  const char *,
  void *
);

//##############################################################################
// main().
//##############################################################################

int
main( int argc, char **argv )
{

  Libnucnet *p_my_nucnet;

  if( argc != 3 )
  {
    fprintf(
      stderr,
      "\nUsage: %s xml_file zone_xpath nuc_xpath\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  xml_file = network xml output file\n\n"
    );
    fprintf(
      stderr,
      "  zone_xpath = XPath expression to select zones\n\n"
    );
    return EXIT_FAILURE;
  }

  p_my_nucnet =
    Libnucnet__new_from_xml(
      argv[1],
      NULL,
      NULL,
      argv[2]
    );

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function)
      nnt::zone_compare_by_first_label
  );

  BOOST_FOREACH( nnt::Zone zone, nnt::make_zone_list( p_my_nucnet ) )
  {

    fprintf(
      stdout,
      "Zone: %s %s %s\n",
      Libnucnet__Zone__getLabel( zone.getNucnetZone(), 1 ),
      Libnucnet__Zone__getLabel( zone.getNucnetZone(), 2 ),
      Libnucnet__Zone__getLabel( zone.getNucnetZone(), 3 )
    );

    Libnucnet__Zone__iterateOptionalProperties(
      zone.getNucnetZone(),
      NULL,
      NULL,
      NULL,
      (Libnucnet__Zone__optional_property_iterate_function)
         print_optional_properties,
      NULL
    );

    std::cout << std::endl;

  }

  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}

//##############################################################################
// print_optional_properties().
//##############################################################################

void
print_optional_properties(
  const char * s_name,
  const char * s_tag1,
  const char * s_tag2,
  const char * s_value,
  void * p_data
)
{

  if( p_data )
  {
    std::cerr
      << "No extra data should be supplied to this routine." << std::endl;
    exit( EXIT_FAILURE );
  }

  if( s_tag1 )
  {
    if( s_tag2 )
    {
      fprintf(
        stdout,
        "Name: %30s  Tag 1: %5s  Tag 2: %5s  Value: %s\n",
        s_name,
        s_tag1,
        s_tag2,
        s_value
      );
    }
    else
    {
      fprintf(
        stdout,
        "Name: %30s  Tag 1: %5s  Value: %s\n",
        s_name,
        s_tag1,
        s_value
      );
    }
  }
  else
  {
    fprintf(
      stdout,
      "Name: %30s  Value: %s\n",
      s_name,
      s_value
    );
  }

}
