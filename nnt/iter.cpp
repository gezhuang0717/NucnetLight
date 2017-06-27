////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2012 Clemson University.
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
//! \brief Code for iteration of zones, reactions, species, reactants, and
//!    products.
//!
////////////////////////////////////////////////////////////////////////////////

#include "nnt/iter.h"

/**
 * @brief The NucNet Tools namespace.
 */
namespace nnt
{

//##############################################################################
// insert_zone_in_list().
//##############################################################################

int
insert_zone_in_list(
  Libnucnet__Zone * p_zone,
  zone_list_t * p_list
)
{

  Zone * pNewZone = new Zone;

  pNewZone->setNucnetZone( p_zone );

  p_list->push_back( pNewZone );

  return 1;

}

//##############################################################################
// insert_species_in_list().
//##############################################################################

int
insert_species_in_list(
  Libnucnet__Species * p_species,
  species_list_t * p_list
)
{

  Species * pNewSpecies = new Species;

  pNewSpecies->setNucnetSpecies( p_species );

  p_list->push_back( pNewSpecies );

  return 1;

}

//##############################################################################
// insert_reaction_in_list().
//##############################################################################

int
insert_reaction_in_list(
  Libnucnet__Reaction * p_reaction,
  reaction_list_t * p_list
)
{

  Reaction * pNewReaction = new Reaction;

  pNewReaction->setNucnetReaction( p_reaction );

  p_list->push_back( pNewReaction );

  return 1;

}

//##############################################################################
// insert_reaction_element_in_list().
//##############################################################################

int
insert_reaction_element_in_list(
  Libnucnet__Reaction__Element * p_reaction_element,
  reaction_element_list_t * p_list
)
{

  ReactionElement * pNewReactionElement = new ReactionElement;

  pNewReactionElement->setNucnetReactionElement( p_reaction_element );

  p_list->push_back( pNewReactionElement );

  return 1;

}

//##############################################################################
// make_zone_list().
//##############################################################################

/**
 * \brief Make a list of Zones from the input Libnucnet structure.
 * \param p_my_nucnet A pointer to a Libnucnet structure.
 * \return A pointer to a new boost::ptr_list of Zones.
 */

zone_list_t make_zone_list( Libnucnet * p_my_nucnet )
{

  zone_list_t my_list;

  Libnucnet__iterateZones(
    p_my_nucnet,
    (Libnucnet__Zone__iterateFunction) insert_zone_in_list,
    &my_list
  );

  return my_list;

}

//##############################################################################
// make_species_list().
//##############################################################################

/**
 * \brief Make a list of Species from the input Libnucnet__Nuc structure.
 * \param p_my_nuc A pointer to a Libnucnet__Nuc structure.
 * \return The populated boost::ptr_list of Species.
 */

species_list_t make_species_list( Libnucnet__Nuc * p_my_nuc )
{

  species_list_t my_list;

  Libnucnet__Nuc__iterateSpecies(
    p_my_nuc,
    (Libnucnet__Species__iterateFunction) insert_species_in_list,
    &my_list
  );

  return my_list;

}

//##############################################################################
// make_reaction_list().
//##############################################################################

/**
 * \brief Make a list of Reactions from the input Libnucnet__Reac structure.
 * \param p_my_reac A pointer to a Libnucnet__Reac structure.
 * \return A pointer to a new boost::ptr_list of Reactions.
 */

reaction_list_t make_reaction_list( Libnucnet__Reac * p_my_reac )
{

  reaction_list_t my_list;

  Libnucnet__Reac__iterateReactions(
    p_my_reac,
    (Libnucnet__Reaction__iterateFunction) insert_reaction_in_list,
    &my_list
  );

  return my_list;

}

//##############################################################################
// make_reaction_reactant_list().
//##############################################################################

/**
 * \brief Make a list of all reactants in a reaction from the input
	  Libnucnet__Reaction structure.
 * \param p_my_reaction A pointer to a Libnucnet__Reaction structure.
 * \return A pointer to a new boost::ptr_list of the reactants in a reaction.
 */

reaction_element_list_t
make_reaction_reactant_list( Libnucnet__Reaction * p_my_reaction )
{

  reaction_element_list_t my_list;

  Libnucnet__Reaction__iterateReactants(
    p_my_reaction,
    (Libnucnet__Reaction__Element__iterateFunction)
      insert_reaction_element_in_list,
    &my_list
  );

  return my_list;

}

//##############################################################################
// make_reaction_nuclide_reactant_list().
//##############################################################################

/**
 * \brief Make a list of nuclide reactants in a reaction from the input
	  Libnucnet__Reaction structure.
 * \param p_my_reaction A pointer to a Libnucnet__Reaction structure.
 * \return A pointer to a new boost::ptr_list of the reactants in a reaction.
 */

reaction_element_list_t
make_reaction_nuclide_reactant_list( Libnucnet__Reaction * p_my_reaction )
{

  reaction_element_list_t my_list;

  Libnucnet__Reaction__iterateNuclideReactants(
    p_my_reaction,
    (Libnucnet__Reaction__Element__iterateFunction)
      insert_reaction_element_in_list,
    &my_list
  );

  return my_list;

}

//##############################################################################
// make_reaction_product_list().
//##############################################################################

/**
 * \brief Make a list of products in a reaction from the input
	  Libnucnet__Reaction structure.
 * \param p_my_reaction A pointer to a Libnucnet__Reaction structure.
 * \return A pointer to a new boost::ptr_list of the products in a reaction.
 */

reaction_element_list_t
make_reaction_product_list( Libnucnet__Reaction * p_my_reaction )
{

  reaction_element_list_t my_list;

  Libnucnet__Reaction__iterateProducts(
    p_my_reaction,
    (Libnucnet__Reaction__Element__iterateFunction)
      insert_reaction_element_in_list,
    &my_list
  );

  return my_list;

}

//##############################################################################
// make_reaction_nuclide_product_list().
//##############################################################################

/**
 * \brief Make a list of nuclide products in a reaction from the input
	  Libnucnet__Reaction structure.
 * \param p_my_reaction A pointer to a Libnucnet__Reaction structure.
 * \return A pointer to a new boost::ptr_list of the products in a reaction.
 */

reaction_element_list_t
make_reaction_nuclide_product_list( Libnucnet__Reaction * p_my_reaction )
{

  reaction_element_list_t my_list;

  Libnucnet__Reaction__iterateNuclideProducts(
    p_my_reaction,
    (Libnucnet__Reaction__Element__iterateFunction)
      insert_reaction_element_in_list,
    &my_list
  );

  return my_list;

}

} // namespace nnt
