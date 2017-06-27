////////////////////////////////////////////////////////////////////////////////
// This file was originally written by Tianhong Yu and Bradley S. Meyer.
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
//! \brief Example code to compute the flows into and out of a cluster.
////////////////////////////////////////////////////////////////////////////////

#include <boost/format.hpp>
#include <boost/tuple/tuple.hpp>

#include <Libnucnet.h>
#include "nnt/auxiliary.h"
#include "nnt/iter.h"
#include "nnt/string_defs.h"
#include "nnt/two_d_weak_rates.h"
#include "nnt/weak_detailed_balance.h"

#include "user/user_rate_functions.h"
#include "user/flow_utilities.h"
#include "user/network_utilities.h"

int main( int argc, char * argv[] ) {

  Libnucnet *p_my_nucnet;
  Libnucnet__NetView * p_net_view, * p_cluster_view;
  std::pair<double,double> flows;
  double d_min;

  //============================================================================
  // Check input.
  //============================================================================

   if ( argc < 4 || argc > 6 ) {
      fprintf(
        stderr,
        "\nUsage: %s xml_filename zone_xpath cluster_xpath min_flow reac_xpath\n\n",
        argv[0]
      );
      fprintf(
        stderr, "  xml_filename = input xml filename\n\n"
      );
      fprintf(
        stderr, "  zone_xpath = XPath to select zones for flows\n\n"
      );
      fprintf(
        stderr, "  cluster_xpath = XPath to select cluster\n\n"
      );
      fprintf(
        stderr,
        "  min_flow = minimum net flow for printout (optional--default=0)\n\n"
      );
      fprintf(
        stderr, "  reac_xpath = reaction xpath (optional)\n\n"
      );
      return EXIT_FAILURE;
   }

  //============================================================================
  // Read file and exit if not present.
  //============================================================================

  p_my_nucnet = Libnucnet__new_from_xml( argv[1], NULL, NULL, argv[2] );

  if( !p_my_nucnet ) {
    fprintf( stderr, "Input data not read!\n" );
    return EXIT_FAILURE;
  }

  //============================================================================
  // Register rate functions.
  //============================================================================

  user::register_my_rate_functions(
    Libnucnet__Net__getReac(
      Libnucnet__getNet( p_my_nucnet )
    )
  );

  //============================================================================
  // Get the minimum net flow for printout.
  //============================================================================

  if( argc > 4 )
    d_min = atof( argv[4] );
  else
    d_min = 0.;

  //============================================================================
  // Get printout formats.
  //============================================================================

  boost::format fmt1(
    "time(s) = %9.4e  t9 = %6.4f  rho(g/cc) = %9.4e  Ye = %6.4f"
  );

  boost::format fmt2( "time(s) = %9.4e  t9 = %6.4f  rho(g/cc) = %9.4e" );

  boost::format fmt3( "%-55s%12.3e%12.3e%12.3e" );

  //============================================================================
  // Get valid net view.
  //============================================================================

  if( argc != 6 )
    p_net_view =
      Libnucnet__NetView__new( Libnucnet__getNet( p_my_nucnet ), "", "" );
  else
    p_net_view =
      Libnucnet__NetView__new( Libnucnet__getNet( p_my_nucnet ), "", argv[5] );

  //============================================================================
  // Set compare function.
  //============================================================================

  Libnucnet__Reac__setReactionCompareFunction(
    Libnucnet__Net__getReac( Libnucnet__NetView__getNet( p_net_view ) ),
    (Libnucnet__Reaction__compare_function)
       nnt::compare_reactions_by_string
  );

  //============================================================================
  // Get the reactions.
  //============================================================================

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list(
      Libnucnet__Net__getReac( Libnucnet__NetView__getNet( p_net_view ) )
    );

  //============================================================================
  // Get the cluster view.
  //============================================================================

  p_cluster_view =
    Libnucnet__NetView__new(
      Libnucnet__getNet( p_my_nucnet ),
      argv[3],
      ""
    );

  //============================================================================
  // Iterate the zones.
  //============================================================================

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  BOOST_FOREACH( nnt::Zone zone, nnt::make_zone_list( p_my_nucnet ) )
  {

    Libnucnet__Zone__updateNetView(
      zone.getNucnetZone(),
      argv[3],
      "",
      NULL,
      Libnucnet__NetView__copy( p_cluster_view )
    );
      
    user::update_my_rate_functions_data( zone );

    //==========================================================================
    // Print conditions.
    //==========================================================================

    if( zone.hasProperty( nnt::s_TIME ) )
    {
      fmt1 %
        boost::lexical_cast<double>( zone.getProperty( nnt::s_TIME ) ) %
        boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) %
        boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) %
        Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 );
      std::cout << fmt1.str() << std::endl;
    }
    else
    {
      fmt1 %
        boost::lexical_cast<double>( zone.getProperty( nnt::s_TIME ) ) %
        boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) %
        boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) );
      std::cout << fmt2.str() << std::endl;
    }

    //==========================================================================
    // Print flows.
    //==========================================================================
  
    double d_forward_in = 0.;
    double d_forward_out = 0.;
    double d_reverse_in = 0.;
    double d_reverse_out = 0.;

    std::vector< boost::tuple<std::string,double,double> > reactions;

    BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
    {

      int i_nuc = 0;

      nnt::reaction_element_list_t reactant_list =
        nnt::make_reaction_nuclide_reactant_list(
          reaction.getNucnetReaction()
        );

      BOOST_FOREACH( nnt::ReactionElement reactant, reactant_list )
      {
        if(
          Libnucnet__Nuc__getSpeciesByName(
            Libnucnet__Net__getNuc(
              Libnucnet__NetView__getNet( p_cluster_view )
            ),
            Libnucnet__Reaction__Element__getName(
              reactant.getNucnetReactionElement()
            )
          )
        )
          i_nuc--;
      }

      nnt::reaction_element_list_t product_list =
        nnt::make_reaction_nuclide_product_list(
          reaction.getNucnetReaction()
        );

      BOOST_FOREACH( nnt::ReactionElement product, product_list )
      {
        if(
          Libnucnet__Nuc__getSpeciesByName(
            Libnucnet__Net__getNuc(
              Libnucnet__NetView__getNet( p_cluster_view )
            ),
            Libnucnet__Reaction__Element__getName(
              product.getNucnetReactionElement()
            )
          )
        )
          i_nuc++;
      }

      if( i_nuc != 0 )
      {

        flows =
          user::compute_flows_for_reaction(
            zone,
            reaction.getNucnetReaction()
          );

        reactions.push_back(
          boost::make_tuple(
            Libnucnet__Reaction__getString( reaction.getNucnetReaction() ),
            flows.first,
            flows.second
          )
        );

        if( i_nuc > 0 )
        {
          d_forward_in += flows.first;
          d_reverse_out += flows.second;
        }
        else
        {
          d_forward_out += flows.first;
          d_reverse_in += flows.second;
        }

      }

    }

    std::cout << "Ycdot_forward_in = " << d_forward_in << std::endl;
    std::cout << "Ycdot_forward_out = " << d_forward_out << std::endl;
    std::cout << "Ycdot_reverse_in = " << d_reverse_in << std::endl;
    std::cout << "Ycdot_reverse_out = " << d_reverse_out << std::endl;
    std::cout <<
      "Ycdot_total_in = " << d_forward_in + d_reverse_in << std::endl;
    std::cout <<
      "Ycdot_total_out = " << d_forward_out + d_reverse_out << std::endl;
    std::cout << "Ycdot = " <<
      d_forward_in - d_forward_out + d_reverse_in - d_reverse_out <<
      std::endl;

    std::cout << "Yc = " <<
      user::compute_cluster_abundance_moment( zone, argv[3], "z", 0  ) <<
      std::endl;

    fprintf( stdout, "\n\t\t\tReaction\t\t\t  Forward     Reverse     Net\n" );
    printf(
      "=======================================================   =========   =========   =========\n"
    );

    for( size_t i = 0; i < reactions.size(); i++ )
    {
      if( fabs( reactions[i].get<1>() - reactions[i].get<2>() ) >= d_min )
      {
        fmt3 %
	  reactions[i].get<0>().c_str() %
	  reactions[i].get<1>() %
	  reactions[i].get<2>() %
	  (reactions[i].get<1>() - reactions[i].get<2>());
        std::cout << fmt3.str() << std::endl;
      }
    }

    std::cout << std::endl;

  }

  //===========================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__NetView__free( p_net_view );
  Libnucnet__NetView__free( p_cluster_view );
  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
