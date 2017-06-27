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
//! \brief Helper header file for scaling routine.
////////////////////////////////////////////////////////////////////////////////

#ifndef SCALING_H
#define SCALING_H

//##############################################################################
// Includes.
//##############################################################################

#include <vector>
#include "nnt/graph.h"

#include "graph_types.h"

#include "graph_helper.h"

//##############################################################################
// Prototypes.
//##############################################################################

void
scale_graph_weights(
  my_graph_t& g,
  std::string
);

#endif  // SCALING_H
