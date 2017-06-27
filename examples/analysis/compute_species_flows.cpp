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
//! \brief Example code to compute flows and steady-state abundance for
//!    zones selected by XPath.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <string>
#include <iostream>
#include <Libnucnet.h>

#include <boost/format.hpp>

#include "nnt/auxiliary.h"
#include "nnt/iter.h"

#include "user/flow_utilities.h"

#include "user/user_rate_functions.h"

#ifdef MY_USER
#include "my_user/my_rates.h"
#endif

//##############################################################################
// check_input().
//##############################################################################

void
check_input( int argc, char **argv )
{

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {
    std::cout << std::endl;
    std::cout << argv[0] <<
      " my_output.xml n15 \"[position() >= last() - 5]\"" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc < 4 || argc > 5 )
  {

    std::cout << std::endl;
    std::cout << "Purpose: " << argv[0] <<
      " prints out production and destruction flows" <<
      std::endl <<
      " for a species from a network xml file for selected zones and selected"
      << std::endl <<
      " reactions." <<
      std::endl;
    fprintf(
      stderr,
      "\nUsage: %s xml_file species_name zone_xpath reac_xpath\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  xml_file = network xml output file\n\n"
    );
    fprintf(
      stderr,
      "  species_name = species for steady flow\n\n"
    );
    fprintf(
      stderr,
      "  zone_xpath = XPath expression to select zone\n\n"
    );
    fprintf(
      stderr,
      "  reac_xpath = XPath expression to select reactions (optional)\n\n"
    );
    std::cout << "For an example usage, type " << std::endl << std::endl;
    std::cout << argv[0] << " --example" << std::endl << std::endl;
    exit( EXIT_FAILURE );

  }

}

//##############################################################################
// add_mod_views().  This copies the views from the old zone to the new.
//##############################################################################

void
add_mod_views(
  nnt::Zone& old_zone,
  nnt::Zone& new_zone
)
{

  user::view_multi views;

  Libnucnet__Zone__iterateOptionalProperties(
    old_zone.getNucnetZone(),
    NULL,
    nnt::s_RATE_MODIFICATION_VIEW,
    NULL,
    (Libnucnet__Zone__optional_property_iterate_function)
      user::modify_rates_views,
    &views
  );

  BOOST_FOREACH( user::view my_view, views )
  {

    Libnucnet__Zone__updateNetView(
      new_zone.getNucnetZone(),
      my_view.nuc_xpath.c_str(),
      my_view.reac_xpath.c_str(),
      NULL,
      Libnucnet__NetView__copy(
        old_zone.getNetView(
          my_view.nuc_xpath.c_str(),
          my_view.reac_xpath.c_str()
        )
      )
    );

  }

}
          
//##############################################################################
// main().
//##############################################################################

int
main( int argc, char **argv )
{

  Libnucnet *p_my_nucnet;
  Libnucnet__ReacView * p_destroy, * p_produce;
  Libnucnet__Species * p_species;
  std::string s_destroy, s_produce;

  check_input( argc, argv );

  if( argc == 4 )
    p_my_nucnet =
      Libnucnet__new_from_xml(
        argv[1],
        NULL,
        NULL,
        argv[3]
      );
  else
    p_my_nucnet =
      Libnucnet__new_from_xml(
        argv[1],
        NULL,
        argv[4],
        argv[3]
      );

  //============================================================================
  // Register user-supplied rate functions and set data.
  //============================================================================
  
  user::register_my_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );
  
#ifdef MY_USER
  my_user::register_my_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );
#endif

  //============================================================================
  // Select species.
  //============================================================================
  
  p_species =
    Libnucnet__Nuc__getSpeciesByName(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
      argv[2]
    );

  //============================================================================
  // Get destruction view.
  //============================================================================
  
  s_destroy =
    std::string( "[reactant = '" ) +
    std::string( argv[2] ) +
    std::string( "']" );

  p_destroy =
    Libnucnet__ReacView__new(
      Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) ),
      s_destroy.c_str()
    );

  nnt::reaction_list_t list1 =
    nnt::make_reaction_list(
      Libnucnet__ReacView__getReac( p_destroy )
    );

  //============================================================================
  // Get production view.
  //============================================================================
  
  s_produce =
    std::string( "[product = '" ) +
    std::string( argv[2] ) +
    std::string( "']" );

  p_produce =
    Libnucnet__ReacView__new(
      Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) ),
      s_produce.c_str()
    );

  nnt::reaction_list_t list2 =
    nnt::make_reaction_list(
      Libnucnet__ReacView__getReac( p_produce )
    );

  //============================================================================
  // Iterate reactions.
  //============================================================================
  
  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function)
      nnt::zone_compare_by_first_label
  );

  nnt::zone_list_t zone_list = nnt::make_zone_list( p_my_nucnet );

  nnt::Zone first_zone = *(zone_list.begin());

  BOOST_FOREACH( nnt::Zone zone, zone_list )
  {

    add_mod_views( first_zone, zone ); 

    double d_destroy_1 = 0, d_destroy_2 = 0, d_destroy_net = 0;

    BOOST_FOREACH( nnt::Reaction reaction, list1 )
    {

      std::pair<double, double> flows =
        user::compute_flows_for_reaction(
          zone,
          reaction.getNucnetReaction()
        );

      d_destroy_1 += flows.first;

      d_destroy_2 += flows.second;

      d_destroy_net += flows.first - flows.second;

    }

    double d_produce_1 = 0, d_produce_2 = 0, d_produce_net = 0;

    BOOST_FOREACH( nnt::Reaction reaction, list2 )
    {

      std::pair<double, double> flows =
        user::compute_flows_for_reaction(
          zone,
          reaction.getNucnetReaction()
        );

      d_produce_1 += flows.first;

      d_produce_2 += flows.second;

      d_produce_net += flows.first - flows.second;

    }

    std::cout <<
      boost::format(
        "%.4e  %.4e  %.4e  %.4e  %.4e  %.4e  %.4e  %.4e  %.4e\n"
      ) %
      boost::lexical_cast<double>( zone.getProperty( nnt::s_TIME ) ) %
      d_produce_1 %
      d_produce_2 %
      d_produce_net %
      d_destroy_1 %
      d_destroy_2 %
      d_destroy_net %
      (
        (
          ( d_produce_1 + d_destroy_2 ) /
          ( d_destroy_1 + d_produce_2 + 1.e-300)
        ) *
        Libnucnet__Zone__getSpeciesAbundance(
          zone.getNucnetZone(),
          p_species
        ) 
      ) %
      Libnucnet__Zone__getSpeciesAbundance(
        zone.getNucnetZone(),
        p_species
      );

  }

  Libnucnet__ReacView__free( p_destroy );
  Libnucnet__ReacView__free( p_produce );
  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}


