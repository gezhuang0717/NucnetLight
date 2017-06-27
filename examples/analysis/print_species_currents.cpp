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
//! \brief Example code to write out the integrated currents into and out of
//!        the chosen species.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Include.
//##############################################################################

#include <iostream>

#include <boost/multi_index_container.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/identity.hpp>
#include <boost/multi_index/member.hpp>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

#include "nnt/auxiliary.h"
#include "nnt/iter.h"
#include "nnt/string_defs.h"

//##############################################################################
// Hash to handle reaction currents.
//##############################################################################

struct Reaction_Current
{

  std::string name;
  double current;

  Reaction_Current( const std::string& s, double d ) :
    name( s ), current( d ) {}

  bool operator<(const Reaction_Current &p)
    const { return current > p.current; }

};

typedef boost::multi_index::multi_index_container<
  Reaction_Current,
  boost::multi_index::indexed_by<
    boost::multi_index::ordered_unique<
      boost::multi_index::identity<Reaction_Current>
    >,
    boost::multi_index::hashed_unique<
      boost::multi_index::member<
        Reaction_Current, std::string, &Reaction_Current::name
      >
    >
  >
> current_multi;

//##############################################################################
// main().
//##############################################################################

int
main( int argc, char **argv )
{

  Libnucnet *p_my_nucnet;
  Libnucnet__ReacView * p_view;
  nnt::Zone flow_current_zone;
  std::string s_reac_xpath;
  double d_produce = 0, d_destroy = 0;
  current_multi producing_currents, destroying_currents;

  if( argc != 3 )
  {
    fprintf(
      stderr,
      "\nUsage: %s xml_file species\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  xml_file = network xml file\n\n"
    );
    fprintf(
      stderr,
      "  species = name of the species\n\n"
    );
    return EXIT_FAILURE;
  }

  p_my_nucnet =
    Libnucnet__new_from_xml(
      argv[1],
      NULL,
      NULL,
      NULL
    );

  //============================================================================
  // Get currents zone.
  //============================================================================
  
  Libnucnet__Zone * p_currents_zone =
    Libnucnet__getZoneByLabels( p_my_nucnet, "integrated currents", "0", "0" );

  if( !p_currents_zone )
  {
    std::cerr << "No integrated currents zone in " << argv[1] << std::endl;
  }

  flow_current_zone.setNucnetZone( p_currents_zone );

  //============================================================================
  // Production reactions.
  //============================================================================
  
  std::cout <<
    boost::format( "\n\t\t\tProduction Reaction\t\t Current\n" );

  std::cout <<
    boost::format(
      "=======================================================\t ==========\n"
    );

  s_reac_xpath = std::string( "[product = '" ) + argv[2] + std::string( "']" );

  p_view =
    Libnucnet__ReacView__new(
      Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) ),
      s_reac_xpath.c_str()
    );

  nnt::reaction_list_t produce_list =
    nnt::make_reaction_list( Libnucnet__ReacView__getReac( p_view ) );

  BOOST_FOREACH( nnt::Reaction reaction, produce_list )
  {

    double d_current =
      static_cast<double>(
        nnt::reaction_element_count(
          reaction,
          "product",
          argv[2]
        )
      ) *
      boost::lexical_cast<double>(
        flow_current_zone.getProperty(
          nnt::s_FLOW_CURRENT,
          Libnucnet__Reaction__getString( reaction.getNucnetReaction() )
        )
      );

    producing_currents.insert(
      Reaction_Current(
        Libnucnet__Reaction__getString( reaction.getNucnetReaction() ),
        d_current
      )
    );

    d_produce += d_current;

  }
      
   BOOST_FOREACH( Reaction_Current my_current, producing_currents )
  {

    std::cout <<
      boost::format( "%-55s%12.4e\n" ) %
      my_current.name.c_str() %
      my_current.current;

  }

  Libnucnet__ReacView__free( p_view );

  //============================================================================
  // Destruction reactions.
  //============================================================================
  
  std::cout <<
    boost::format( "\n\t\t\tDestruction Reaction\t\t Current\n" );

  std::cout <<
    boost::format(
      "=======================================================\t ==========\n"
    );

  s_reac_xpath = std::string( "[reactant = '" ) + argv[2] + std::string( "']" );

  p_view =
    Libnucnet__ReacView__new(
      Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) ),
      s_reac_xpath.c_str()
    );

  nnt::reaction_list_t destroy_list =
    nnt::make_reaction_list( Libnucnet__ReacView__getReac( p_view ) );

  BOOST_FOREACH( nnt::Reaction reaction, destroy_list )
  {

    double d_current =
      static_cast<double>(
        nnt::reaction_element_count(
          reaction,
          "reactant",
          argv[2]
        )
      ) *
      boost::lexical_cast<double>(
        flow_current_zone.getProperty(
          nnt::s_FLOW_CURRENT,
          Libnucnet__Reaction__getString( reaction.getNucnetReaction() )
        )
      );

    destroying_currents.insert(
      Reaction_Current(
        Libnucnet__Reaction__getString( reaction.getNucnetReaction() ),
        d_current
      )
    );

    d_destroy += d_current;

  }
      
  Libnucnet__ReacView__free( p_view );

  BOOST_FOREACH( Reaction_Current my_current, destroying_currents )
  {

    std::cout <<
      boost::format( "%-55s%12.4e\n" ) %
      my_current.name.c_str() %
      my_current.current;

  }

  //============================================================================
  // Total diagnostics.
  //============================================================================

  std::cout <<
    boost::format( "\nTotal integrated production current: %12.4e\n" ) %
    d_produce;

  std::cout <<
    boost::format( "Total integrated destruction current: %12.4e\n" ) %
    d_destroy;

  std::cout <<
    boost::format( "Total integrated net current: %12.4e\n" ) %
    (d_produce - d_destroy);

  double d_initial =
    boost::lexical_cast<double>(
      flow_current_zone.getProperty( 
        nnt::s_INITIAL_ABUNDANCE,
        argv[2]
      )
    );

  double d_final =
    boost::lexical_cast<double>(
      flow_current_zone.getProperty( 
        nnt::s_FINAL_ABUNDANCE,
        argv[2]
      )
    );

  std::cout <<
    boost::format( "\nInitial abundance of %s: %12.4e\n" ) %
    argv[2] %
    d_initial;

  std::cout <<
    boost::format( "Final abundance of %s: %12.4e\n" ) %
    argv[2] %
    d_final;

  std::cout <<
    boost::format( "Abundance change of %s: %12.4e\n\n" ) %
    argv[2] %
    (d_final - d_initial);

  //============================================================================
  // Clean up.
  //============================================================================

  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
