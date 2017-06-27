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
//! \brief Example code to compute the flows in given zone (chosen by
//!    XPath) by nucleon number.
////////////////////////////////////////////////////////////////////////////////

#include <Libnucnet.h>

#include <boost/format.hpp>

#include "nnt/auxiliary.h"
#include "nnt/iter.h"
#include "nnt/string_defs.h"
#include "nnt/two_d_weak_rates.h"
#include "nnt/weak_detailed_balance.h"

#include "user/user_rate_functions.h"
#include "user/flow_utilities.h"

#ifdef MY_USER
#include "my_user/my_rates.h"
#endif

std::vector<nnt::Reaction>
get_nucleon_number_changing_reactions(
  Libnucnet__Net * p_net,
  std::string nucleon
)
{

  std::vector<nnt::Reaction> reactions;

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list(
      Libnucnet__Net__getReac( p_net )
    );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    int delta = 0;

    nnt::reaction_element_list_t reactant_list =
      nnt::make_reaction_nuclide_reactant_list(
	reaction.getNucnetReaction()
      );
       
    BOOST_FOREACH( nnt::ReactionElement reactant, reactant_list )
    {

      Libnucnet__Species * p_species =
	Libnucnet__Nuc__getSpeciesByName(
	  Libnucnet__Net__getNuc( p_net ),
	  Libnucnet__Reaction__Element__getName(
	    reactant.getNucnetReactionElement()
	  )
	);

      if( Libnucnet__Species__getZ( p_species ) > 2 )
      {

	if( nucleon == "z" )
	  delta += Libnucnet__Species__getZ( p_species );
	else if( nucleon == "n" )
	  delta += Libnucnet__Species__getA( p_species ) -
		  Libnucnet__Species__getZ( p_species );
	else if( nucleon == "a" )
	  delta += Libnucnet__Species__getA( p_species );

      }

    }

    nnt::reaction_element_list_t product_list =
      nnt::make_reaction_nuclide_product_list(
	reaction.getNucnetReaction()
      );

    BOOST_FOREACH( nnt::ReactionElement product, product_list )
    {

      Libnucnet__Species * p_species =
	Libnucnet__Nuc__getSpeciesByName(
	  Libnucnet__Net__getNuc( p_net ),
	  Libnucnet__Reaction__Element__getName(
	    product.getNucnetReactionElement()
	  )
	);

      if( Libnucnet__Species__getZ( p_species ) > 2 )
      {

	if( nucleon == "z" )
	  delta -= Libnucnet__Species__getZ( p_species );
	else if( nucleon == "n" )
	  delta -= Libnucnet__Species__getA( p_species ) -
		  Libnucnet__Species__getZ( p_species );
	else if( nucleon == "a" )
	  delta -= Libnucnet__Species__getA( p_species );

      }

    }

    if( delta != 0 ) reactions.push_back( reaction );

  }

  return reactions;

}

int main( int argc, char * argv[] ) {

  Libnucnet *p_my_nucnet;
  Libnucnet__NetView * p_net_view;
  std::vector<nnt::Reaction> reactions;
  std::pair<double,double> flows;
  size_t i;

  //============================================================================
  // Check input.
  //============================================================================

   if ( argc < 4 || argc > 5 ) {
      fprintf(
        stderr,
        "\nUsage: %s xml_filename zone_xpath nucleon_type reac_xpath min_flow\n\n",
        argv[0]
      );
      fprintf(
        stderr, "  xml_filename = input xml filename\n\n"
      );
      fprintf(
        stderr, "  zone_xpath = XPath to select zones for flows\n\n"
      );
      fprintf(
        stderr, "  nucleon_type = flag to choose nucleon type (z,n,a)\n\n"
      );
      fprintf(
        stderr, "  reac_xpath = reaction xpath (optional--required if min_flow present)\n\n"
      );
      fprintf(
        stderr, "  min_flow = minimum net flow to print out (optional--if not set, min is 1.e-50)\n\n"
      );
      return EXIT_FAILURE;
   }

  //============================================================================
  // Check nucleon type.
  //============================================================================

  if(
    std::string( argv[3] ) != "z" ||
    std::string( argv[3] ) != "n" ||
    std::string( argv[3] ) != "a"
  )
  {
    std::cerr << "Nucleon type must be z, n, or a." << std::endl;
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

#ifdef MY_USER
  my_user::register_my_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );
#endif

  //============================================================================
  // Get the zones.
  //============================================================================

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  nnt::zone_list_t zone_list = nnt::make_zone_list( p_my_nucnet );

  //============================================================================
  // Get valid net view.
  //============================================================================

  if( argc == 4 )
    p_net_view =
      Libnucnet__NetView__new( Libnucnet__getNet( p_my_nucnet ), "", "" );
  else
    p_net_view =
      Libnucnet__NetView__new( Libnucnet__getNet( p_my_nucnet ), "", argv[4] );

  //==========================================================================
  // Assign arrays and output largest nucleon number.
  //==========================================================================

  std::size_t i_max =
    Libnucnet__Nuc__getLargestNucleonNumber(
      Libnucnet__Net__getNuc(
        Libnucnet__NetView__getNet( p_net_view )
      ),
      argv[3]
    );

  std::vector<double> in_forward_flows( i_max + 1 );

  std::vector<double> in_reverse_flows( i_max + 1 );

  std::vector<double> out_forward_flows( i_max + 1 );

  std::vector<double> out_reverse_flows( i_max + 1 );

  std::cout << i_max << std::endl;

  //============================================================================
  // Get the reactions.
  //============================================================================

  reactions =
    get_nucleon_number_changing_reactions(
      Libnucnet__NetView__getNet( p_net_view ),
      argv[3]
    );

  Libnucnet__NetView__free( p_net_view );

  //============================================================================
  // Iterate the zones.
  //============================================================================

  BOOST_FOREACH( nnt::Zone zone, zone_list )
  {
   
    if( !( zone == *(zone_list.begin()) ) )
    {
      Libnucnet__Zone__copy_net_views(
        zone.getNucnetZone(),
        (*zone_list.begin()).getNucnetZone()
      );
    }

    user::update_my_rate_functions_data( zone );

    //==========================================================================
    // Print conditions.
    //==========================================================================

    if( zone.hasProperty( nnt::s_TIME ) )
      std::cout <<
        boost::format( "time(s) = %12.4e t9 = %12.4e rho(g/cc) = %12.4e\n" ) %
        ( boost::lexical_cast<double>( zone.getProperty( nnt::s_TIME ) ) ) %
        ( boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) ) %
        ( boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) );
    else
      std::cout <<
        boost::format( "t9 = %12.4e rho(g/cc) = %12.4e\n" ) %
        ( boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) ) %
        ( boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) );

    //==========================================================================
    // Zero out arrays.
    //==========================================================================

    std::fill( in_forward_flows.begin(), in_forward_flows.end(), 0. );
    std::fill( in_reverse_flows.begin(), in_reverse_flows.end(), 0. );
    std::fill( out_forward_flows.begin(), out_forward_flows.end(), 0. );
    std::fill( out_reverse_flows.begin(), out_reverse_flows.end(), 0. );

    //==========================================================================
    // Compute and assign flows.
    //==========================================================================

    BOOST_FOREACH( nnt::Reaction reaction, reactions )
    {

      flows =
	user::compute_flows_for_reaction(
	  zone,
	  reaction.getNucnetReaction()
      );

      nnt::reaction_element_list_t reactant_list =
	nnt::make_reaction_nuclide_reactant_list(
	  reaction.getNucnetReaction()
	);
   
      BOOST_FOREACH( nnt::ReactionElement reactant, reactant_list )
      {

	Libnucnet__Species * p_species =
	  Libnucnet__Nuc__getSpeciesByName(
	    Libnucnet__Net__getNuc(
	      Libnucnet__Zone__getNet( zone.getNucnetZone() )
	    ),
	    Libnucnet__Reaction__Element__getName(
	      reactant.getNucnetReactionElement()
	    )
	  );

	if( std::string( argv[3] ) == "z" )
	  i = Libnucnet__Species__getZ( p_species );
	else if( std::string( argv[3] ) == "n" )
	  i = Libnucnet__Species__getA( p_species ) -
		 Libnucnet__Species__getZ( p_species );
	else if( std::string( argv[3] ) == "a" )
	  i = Libnucnet__Species__getA( p_species );

	out_forward_flows[i] += flows.first;
	in_reverse_flows[i] += flows.second;

      }

      nnt::reaction_element_list_t product_list =
	nnt::make_reaction_nuclide_product_list(
	  reaction.getNucnetReaction()
	);

      BOOST_FOREACH( nnt::ReactionElement product, product_list )
      {

	Libnucnet__Species * p_species =
	  Libnucnet__Nuc__getSpeciesByName(
	    Libnucnet__Net__getNuc(
	      Libnucnet__Zone__getNet( zone.getNucnetZone() )
	    ),
	    Libnucnet__Reaction__Element__getName(
	      product.getNucnetReactionElement()
	    )
	  );

	if( std::string( argv[3] ) == "z" )
	  i = Libnucnet__Species__getZ( p_species );
	else if( std::string( argv[3] ) == "n" )
	  i = Libnucnet__Species__getA( p_species ) -
		 Libnucnet__Species__getZ( p_species );
	else if( std::string( argv[3] ) == "a" )
	  i = Libnucnet__Species__getA( p_species );


	out_reverse_flows[i] += flows.second;
	in_forward_flows[i] += flows.first;

      }

    }

    //==========================================================================
    // Print out.
    //==========================================================================

    for( size_t j = 0; j < in_forward_flows.size(); j++ )
    {

      std::cout <<
        boost::format( "%5d %12.4e %12.4e %12.4e %12.4e\n" ) %
        j %
        in_forward_flows[j] %
        in_reverse_flows[j] %
        out_forward_flows[j] %
        out_reverse_flows[j];

    }

    std::cout << std::endl;

  }

  //===========================================================================
  // Clean up and exit.
  //==========================================================================*/

  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
