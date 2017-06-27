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
//! \brief Helper header file for graph routines.
////////////////////////////////////////////////////////////////////////////////

#ifndef GRAPH_HELPER_H
#define GRAPH_HELPER_H

//##############################################################################
// Includes.
//##############################################################################

#include <vector>
#include "nnt/graph.h"

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/mem_fun.hpp>

#include "graph_types.h"

#include "color.h"
#include "linestyle.h"
#include "scaling.h"

//##############################################################################
// Structure to handle duplicate reactant or product links.
//##############################################################################

struct Graph_Link
{
  std::string out_name;
  std::string in_name;

  Graph_Link( const std::string &n1, const std::string &n2 )
    : out_name(n1), in_name(n2) {}

  bool operator==(const Graph_Link &p) const
  {
    return out_name == p.out_name && in_name == p.in_name;
  }
};

//##############################################################################
// Prototypes.
//##############################################################################

std::size_t hash_value(Graph_Link const &);

my_graph_t
create_induced_subgraph(
  my_graph_t&,
  Libnucnet__NucView *,
  bool
);

void
add_reaction_net_flow_edges(
  my_graph_t&,
  my_vertex_hash_t&,
  nnt::Reaction&,
  std::pair<double,double>&
);

void
add_reaction_all_flow_edges(
  my_graph_t&,
  my_vertex_hash_t&,
  nnt::Reaction&,
  std::pair<double,double>&
);

void
add_solar(
  my_graph_t &,
  Libnucnet__NucView *
);

double
get_max_flow(
  my_graph_t&
);

#endif  // GRAPH_HELPER_H
