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
//! \brief Header file for zone links code.
////////////////////////////////////////////////////////////////////////////////

#ifndef ZONE_LINKS_H
#define ZONE_LINKS_H

//##############################################################################
// Includes.
//##############################################################################

#include <iostream>
#include <string>

#include <Libnucnet.h>

#include "nnt/auxiliary.h"
#include "nnt/containers.h"
#include "nnt/iter.h"
#include "nnt/string_defs.h"

/**
 * @brief A namespace for user-defined rate functions.
 */

namespace user
{

//##############################################################################
// Metafunction for zone links.
//##############################################################################

template <typename T>
class zone_links {
  public:
    typedef
    typename nnt::ordered_bimap_with_data<nnt::Zone,nnt::Zone,T>::type type;
};
 
//##############################################################################
// Prototypes.
//##############################################################################

size_t
count_links( nnt::Zone& );
 
void
links_counter(
  const char *,
  const char *,
  const char *,
  const char *,
  size_t *
);

size_t
count_links( const nnt::Zone& );

void
assign_link_zone_tag(
  const char *,
  const char *,
  const char *,
  const char *,
  boost::tuple<std::string,std::string,std::string>&
);

} // namespace user

#endif // ZONE_LINKS_H
