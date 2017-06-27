////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011 Clemson University.
//
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
// You should have received a copy of the GNU General Public License
// along with this software; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
// USA
//
//////////////////////////////////////////////////////////////////////////////*/

////////////////////////////////////////////////////////////////////////////////
//!
//! \file
//! \brief A header file for iteration of zones,
//!        reactions, species, reactants and products, and properties.
//!
////////////////////////////////////////////////////////////////////////////////

#ifndef NNT_ITER_H
#define NNT_ITER_H

#include <boost/ptr_container/ptr_list.hpp>
#include <boost/foreach.hpp>
#include "nnt/wrappers.h"

/**
 * @brief The NucNet Tools namespace.
 */
namespace nnt
{

typedef boost::ptr_list<Species> species_list_t;

typedef boost::ptr_list<Reaction> reaction_list_t;

typedef boost::ptr_list<ReactionElement> reaction_element_list_t;

typedef boost::ptr_list<Zone> zone_list_t;

zone_list_t make_zone_list( Libnucnet * );

reaction_list_t make_reaction_list( Libnucnet__Reac * );

species_list_t make_species_list( Libnucnet__Nuc * );

reaction_element_list_t
make_reaction_reactant_list( Libnucnet__Reaction * );

reaction_element_list_t
make_reaction_nuclide_reactant_list( Libnucnet__Reaction * );

reaction_element_list_t
make_reaction_product_list( Libnucnet__Reaction * );

reaction_element_list_t
make_reaction_nuclide_product_list( Libnucnet__Reaction * );

} // namespace nnt

#endif // NNT_ITER_HPP
