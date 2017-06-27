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
//! \brief Example code to compute the difference between integrated
//!    currents in two currents xml files.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Include.
//##############################################################################

#include <map>

#include <boost/lexical_cast.hpp>

#include "nnt/iter.h"
#include "nnt/string_defs.h"

//##############################################################################
// global structure.
//##############################################################################

typedef struct {
  Libnucnet__Zone * pZone;
  std::string tag;
} my_struct;

//##############################################################################
// copy_abundance().
//##############################################################################

void
copy_abundance(
  const char * s_name,
  const char * s_tag1,
  const char * s_tag2,
  const char * s_value,
  my_struct * p_my_struct
)
{

  Libnucnet__Zone__updateProperty(
    p_my_struct->pZone,
    p_my_struct->tag.c_str(),
    s_tag1,
    s_tag2,
    s_value
  );
    
}

//##############################################################################
// fill_map().
//##############################################################################

void
fill_map(
  const char * s_name,
  const char * s_tag1,
  const char * s_tag2,
  const char * s_value,
  std::map<std::string,double>& map
)
{

  map[s_tag1] = boost::lexical_cast<double>( s_value );
    
}

//##############################################################################
// main().
//##############################################################################

int
main( int argc, char **argv )
{

  Libnucnet * p_my_nucnet;
  Libnucnet__Zone * p_flow_current_zone_1,
                  * p_flow_current_zone_2,
                  * p_new_zone;
  std::map<std::string,double> currents1, currents2;
  my_struct work;

  if( argc != 4 )
  {
    fprintf(
      stderr,
      "\nUsage: %s xml_file1 xml_file2 new_xml\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  xml_file1 = 1st network currents xml file\n\n"
    );
    fprintf(
      stderr,
      "  xml_file2 = 2nd network currents xml file\n\n"
    );
    fprintf(
      stderr,
      "  new_xml = new currents xml file with the difference in currents\n\n"
    );
    return EXIT_FAILURE;
  }

  p_my_nucnet =
    Libnucnet__new_from_xml(
      argv[1],
      NULL,
      NULL,
      NULL
    );

  p_flow_current_zone_1 =
    Libnucnet__getZoneByLabels( p_my_nucnet, "integrated currents", "0", "0" );

  if( !p_flow_current_zone_1 )
  {
    std::cerr << "No integrated currents zone in " << argv[1] << std::endl;
  }

  Libnucnet__relabelZone(
    p_my_nucnet,
    p_flow_current_zone_1,
    "integrated currents 1",
    "0",
    "0"
  );

  Libnucnet__assignZoneDataFromXml( p_my_nucnet, argv[2], NULL );

  p_flow_current_zone_2 =
    Libnucnet__getZoneByLabels( p_my_nucnet, "integrated currents", "0", "0" );

  if( !p_flow_current_zone_2 )
  {
    std::cerr << "No integrated currents zone in " << argv[2] << std::endl;
  }

  Libnucnet__relabelZone(
    p_my_nucnet,
    p_flow_current_zone_2,
    "integrated currents 2",
    "0",
    "0"
  );

  //============================================================================
  // Create output zone and add.
  //============================================================================
  
  p_new_zone =
    Libnucnet__Zone__new(
      Libnucnet__getNet( p_my_nucnet ),
      "integrated currents",
      "0",
      "0"
    );

  Libnucnet__addZone( p_my_nucnet, p_new_zone );

  //============================================================================
  // Fill currents hashes.
  //============================================================================

  Libnucnet__Zone__iterateOptionalProperties(
    p_flow_current_zone_1,
    nnt::s_FLOW_CURRENT,
    NULL,
    NULL,
    (Libnucnet__Zone__optional_property_iterate_function) fill_map,
    &currents1
  );

  Libnucnet__Zone__iterateOptionalProperties(
    p_flow_current_zone_2,
    nnt::s_FLOW_CURRENT,
    NULL,
    NULL,
    (Libnucnet__Zone__optional_property_iterate_function) fill_map,
    &currents2
  );

  //============================================================================
  // Loop on currents.
  //============================================================================

  for(
     std::map<std::string,double>::iterator it = currents1.begin();
     it != currents1.end();
     it++
  )
  {

    std::map<std::string,double>::iterator it2 =
      currents2.find( it->first );

    if( it2 == currents2.end() )
    {
      std::cerr <<
        "Current for reaction " << it->first.c_str() <<
        " is not present in " << argv[2] << std::endl;
    }

    Libnucnet__Zone__updateProperty(
      p_new_zone,
      nnt::s_FLOW_CURRENT,
      it->first.c_str(),
      NULL,
      boost::lexical_cast<std::string>(
        it->second - it2->second
      ).c_str()
    );

  }
    
  //============================================================================
  // Copy abundances.
  //============================================================================

  work.pZone = p_new_zone;
  work.tag = "initial abundance";

  Libnucnet__Zone__iterateOptionalProperties(
    p_flow_current_zone_2,
    "final abundance",
    NULL,
    NULL,
    (Libnucnet__Zone__optional_property_iterate_function)
       copy_abundance,
    &work
  );

  work.tag = "final abundance";

  Libnucnet__Zone__iterateOptionalProperties(
    p_flow_current_zone_1,
    "final abundance",
    NULL,
    NULL,
    (Libnucnet__Zone__optional_property_iterate_function)
       copy_abundance,
    &work
  );

  //============================================================================
  // Remove old zones.
  //============================================================================

  Libnucnet__removeZone( p_my_nucnet, p_flow_current_zone_1 );

  Libnucnet__removeZone( p_my_nucnet, p_flow_current_zone_2 );

  //============================================================================
  // Write to new file.
  //============================================================================

  Libnucnet__writeToXmlFile( p_my_nucnet, argv[3] );

  //============================================================================
  // Clean up.
  //============================================================================

  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
