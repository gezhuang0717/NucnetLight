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
//! \brief Helper header file for color routines.
////////////////////////////////////////////////////////////////////////////////

#ifndef COLOR_H
#define COLOR_H

//##############################################################################
// Includes.
//##############################################################################

#include <sstream>
#include <vector>
#include "nnt/graph.h"

#include <boost/unordered_map.hpp>
#include <boost/format.hpp>

//##############################################################################
// Prototypes.
//##############################################################################

std::pair<std::string,std::string>
solar_color( void );

boost::unordered_map<std::string,std::string>
get_special_vertex_color_map( void );

void
create_reaction_color_map( Libnucnet__Reac * );

std::string
get_reaction_color( Libnucnet__Reaction * );

void
color_reactions( Libnucnet__Reac *, const char *, const char * );

std::string
get_bounding_color( void );

std::string
get_color_from_int( int );

#endif  // COLOR_H
