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
//! \brief File to add line styles to graph.
////////////////////////////////////////////////////////////////////////////////

#include "linestyle.h"

//##############################################################################
// Global sets for reaction types.
//##############################################################################

boost::unordered_map<std::string,std::string> reaction_linestyle_map;

//##############################################################################
// create_reaction_linestyle_map().
//##############################################################################

void
create_reaction_linestyle_map(
  Libnucnet__Reac * p_reac
)
{

  //============================================================================
  // All reactions start solid.
  //============================================================================
  
  set_reactions_linestyle( p_reac, "", "solid" );

  //============================================================================
  // Here we set line styles for reactions.  Change or add as desired.
  //============================================================================

  //----------------------------------------------------------------------------
  // Weak reactions.
  //----------------------------------------------------------------------------
  
  set_reactions_linestyle(
    p_reac,
    nnt::s_WEAK_XPATH,
    "solid"
  );

}

//##############################################################################
// get_reaction_linestyle().
//##############################################################################

std::string
get_reaction_linestyle( Libnucnet__Reaction * p_reaction )
{

  boost::unordered_map<std::string,std::string>::iterator it;

  it =
    reaction_linestyle_map.find( Libnucnet__Reaction__getString( p_reaction ) );

  if( it == reaction_linestyle_map.end() )
  {
    std::cerr << "Reaction linestyle not found." << std::endl;
    exit( EXIT_FAILURE );
  }

  return it->second;

}

//##############################################################################
// set_reactions_linestyle().
//##############################################################################

void
set_reactions_linestyle(
  Libnucnet__Reac * p_reac,
  const char * s_reac_xpath,
  const char * s_linestyle
)
{

  Libnucnet__ReacView * p_view =
    Libnucnet__ReacView__new(
      p_reac,
      s_reac_xpath
    );

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list( Libnucnet__ReacView__getReac( p_view ) );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {
    reaction_linestyle_map[
      Libnucnet__Reaction__getString( reaction.getNucnetReaction() )
    ] = s_linestyle;
  }

  Libnucnet__ReacView__free( p_view );

} 
