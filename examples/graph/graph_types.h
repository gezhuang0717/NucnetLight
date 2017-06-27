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
//! \brief Header file for useful graph types.
////////////////////////////////////////////////////////////////////////////////

#ifndef GRAPH_TYPES_H
#define GRAPH_TYPES_H

//##############################################################################
// Includes.
//##############################################################################

#include "nnt/graph.h"

//##############################################################################
// Types.
//##############################################################################

struct my_vertex_data
{
  std::string color;
  std::string shape;

  my_vertex_data(
    std::string c = std::string(),
    std::string s = std::string()
  ) : color( c ), shape( s ) {}

};

typedef my_vertex_data my_vertex_t;

typedef std::string my_edge_t;

typedef nnt::net_graph<my_vertex_t,my_edge_t>::type my_graph_t;

typedef nnt::vertex_multi<my_vertex_t,my_edge_t>::type my_vertex_hash_t;

typedef nnt::Graph_Vertex<my_vertex_t,my_edge_t>::type my_graph_vertex_t;

#endif // GRAPH_TYPES_H
