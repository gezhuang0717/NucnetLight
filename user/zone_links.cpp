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
//! \brief Utilities code for zone links.
////////////////////////////////////////////////////////////////////////////////

#include "zone_links.h"

/**
 * @brief A namespace for user-defined rate functions.
 */

namespace user
{

//##############################################################################
// links_counter().
//##############################################################################

void
links_counter(
  const char * s_name,
  const char * s_tag1,
  const char * s_tag2,
  const char * s_value,
  size_t * i_links
)
{

  if( !s_name || !s_tag1 || !s_tag2 || !s_value )
  {
    std::cerr << "Invalid input to links counter." << std::endl;
    exit( EXIT_FAILURE );
  }

  (*i_links)++;

}

//##############################################################################
// count_links().
//##############################################################################

size_t
count_links( nnt::Zone& zone )
{

  size_t i_links = 0;

  Libnucnet__Zone__iterateOptionalProperties(
    zone.getNucnetZone(),
    nnt::s_LINK,
    NULL,
    "0",
    (Libnucnet__Zone__optional_property_iterate_function) links_counter,
    &i_links
  );

  return i_links;

}

//##############################################################################
// assign_link_zone_tag().
//##############################################################################

void
assign_link_zone_tag(
  const char * s_name,
  const char * s_tag1,
  const char * s_tag2,
  const char * s_value,
  boost::tuple<std::string,std::string,std::string>& t
)
{

  if( strcmp( s_tag2, "0" ) == 0 )
    t.get<0>() = s_value;
  else if( strcmp( s_tag2, "1" ) == 0 )
    t.get<1>() = s_value;
  else if( strcmp( s_tag2, "2" ) == 0 )
    t.get<2>() = s_value;
  else
  {
     std::cout << s_tag2 << std::endl;
     std::cerr << "Invalid zone tag." << std::endl;
     exit( EXIT_FAILURE );
  }

}

} // namespace user
