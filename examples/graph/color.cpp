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
//! \brief File to add colors to graph.
////////////////////////////////////////////////////////////////////////////////

#include "color.h"

//##############################################################################
// Global sets for reaction types.
//##############################################################################

boost::unordered_map<std::string,std::string> reaction_color_map;

//##############################################################################
// solar_color().
//##############################################################################

std::pair<std::string,std::string>
solar_color( void )
{

  //============================================================================
  // Set the color of the non-stable species.  Change this as desired.
  //============================================================================
  
  std::string s_non_solar_color = "white";

  //============================================================================
  // Set the color of the stable species.  Change this as desired.
  //============================================================================
  
  std::string s_solar_color = "yellow";

  //============================================================================
  // Return the pair.
  //============================================================================
  
  return std::make_pair( s_non_solar_color, s_solar_color );

}
  
//##############################################################################
// get_special_vertex_color_map().
//##############################################################################

boost::unordered_map<std::string,std::string>
get_special_vertex_color_map( void )
{

  boost::unordered_map<std::string,std::string> my_map;

  //============================================================================
  // Change or add to these as desired.
  //============================================================================
  
  my_map["fe56"] = "yellow";
  my_map["ni60"] = "yellow";

  //============================================================================
  // Return the map.
  //============================================================================
  
  return my_map;

}
  
//##############################################################################
// create_reaction_color_map().
//##############################################################################

void
create_reaction_color_map(
  Libnucnet__Reac * p_reac
)
{

  //============================================================================
  // All reactions start black.
  //============================================================================
  
  color_reactions( p_reac, "", "black" );

  //============================================================================
  // Here we set colors for reactions.  Change or add as desired.
  //============================================================================

  //----------------------------------------------------------------------------
  // Color alpha decays blue.
  //----------------------------------------------------------------------------
  
  color_reactions(
    p_reac,
    "[ count( reactant ) = 1 and product = 'he4' ]",
    "blue"
  );

  //----------------------------------------------------------------------------
  // Color weak reactions red.
  //----------------------------------------------------------------------------
  
  color_reactions(
    p_reac,
    nnt::s_WEAK_XPATH,
    "red"
  );

}

//##############################################################################
// get_reaction_color().
//##############################################################################

std::string
get_reaction_color( Libnucnet__Reaction * p_reaction )
{

  boost::unordered_map<std::string,std::string>::iterator it;

  it =
    reaction_color_map.find( Libnucnet__Reaction__getString( p_reaction ) );

  if( it == reaction_color_map.end() )
  {
    std::cerr << "Reaction color not found." << std::endl;
    exit( EXIT_FAILURE );
  }

  return it->second;

}

//##############################################################################
// color_reactions().
//##############################################################################

void
color_reactions(
  Libnucnet__Reac * p_reac,
  const char * s_reac_xpath,
  const char * s_color
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
    reaction_color_map[
      Libnucnet__Reaction__getString( reaction.getNucnetReaction() )
    ] = s_color;
  }

  Libnucnet__ReacView__free( p_view );

} 

//##############################################################################
// get_bounding_color().
//##############################################################################

std::string
get_bounding_color( void )
{

  //============================================================================
  // Return the color for the bounding box of vertices.
  //============================================================================

  return "blue";

}

//##############################################################################
// get_color_from_int().
//##############################################################################

std::string
get_color_from_int( int i )
{

  std::stringstream s_color;

  if( i >= -255 && i <= 0 )    // Pivot light purple to red.
  {
    s_color.str("");
    s_color << "#FF00";
    s_color << boost::format( "%02x" ) % (255 + i);
  }
  else if( i >= -511 && i < -255 )  // Pivot red to yellow.
  {
    s_color.str("");
    s_color << "#FF";
    s_color << boost::format( "%02x" ) % (-i - 255);
    s_color << "00";
  }
  else if( i >= -767 && i < -511 )  // Pivot yellow to green.
  {
    s_color.str("");
    s_color << "#";
    s_color << boost::format( "%02x" ) % (i + 767);
    s_color << "FF00";
  }
  else if( i <= 255 && i > 0 )  // Pivot light purple to blue.
  {
    s_color.str("");
    s_color << "#";
    s_color << boost::format( "%02x" ) % (255 - i);
    s_color << "00";
    s_color << "FF";
  }
  else if( i < 511 && i > 255 )  // Pivot blue to cyan.
  {
    s_color.str("");
    s_color << "#00";
    s_color << boost::format( "%02x" ) % (i - 255);
    s_color << "FF";
  }
  else
  {
    s_color.str( "#808080" );  // Outside range = grey.
  }

  return s_color.str();

}
