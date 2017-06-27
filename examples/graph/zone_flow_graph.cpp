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
//! \brief Example code to write out the net flow graph of a zone in dot
//!    format.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Include.
//##############################################################################

#include <fstream>
#include <ostream>
#include <iostream>
#include <sstream>
#include <iomanip>

#include <boost/graph/graphviz.hpp>

#include "nnt/iter.h"

#include "user/user_rate_functions.h"

#ifdef MY_USER
#include "my_user/my_rates.h"
#endif

#include "graph_types.h"
#include "graph_helper.h"

#define D_SCALE  100
#define D_MAX_FLOW_SCALE 10

#define B_ALLOW_ISOLATED_VERTICES false

//##############################################################################
// The vertex writer.
//##############################################################################

struct
my_vertex_writer
{
  my_vertex_writer( my_graph_t& g_ ) :           g( g_ ) {};
  template <class Vertex>
  void operator()( std::ostream& out, Vertex v )
  {
    int i_z = Libnucnet__Species__getZ( g[v].getNucnetSpecies() );
    int i_n = Libnucnet__Species__getA( g[v].getNucnetSpecies() ) - i_z;
    char * s_latex_string =
      Libnucnet__Species__createLatexString( g[v].getNucnetSpecies() );
    out <<
      boost::format(
        " [texlbl=\"\\huge{$%s$}\" \
          pos=\"%d,%d!\", \
          style=filled, fillcolor=\"%s\" \
         ]"
      ) % 
        s_latex_string %
        ( D_SCALE * i_n ) %
        ( D_SCALE * i_z ) %
        g[v].getExtraData().color
        << std::endl;
     free( s_latex_string );
  }
  my_graph_t& g;
};

//##############################################################################
// The edge writer.
//##############################################################################

struct
my_edge_writer
{
  my_edge_writer( my_graph_t& g_ ) :
    g( g_ ) {};
  template <class Edge>
  void operator()( std::ostream& out, Edge e )
  {
    out.precision(5);
    std::string s_color = get_reaction_color( g[e].getNucnetReaction() );
    std::string s_style = get_reaction_linestyle( g[e].getNucnetReaction() );
    out <<
      boost::format(
        "[style=\"line width = %gpt, %s\" color = \"%s\"]"
      ) %
        ( D_MAX_FLOW_SCALE * g[e].getWeight() ) %
        s_style %
        s_color
       <<
      std::endl;
  }
  my_graph_t &g;
};

//##############################################################################
// The graph writer.
//##############################################################################

struct
my_graph_writer
{
  my_graph_writer( nnt::Zone& zone_ ) : zone( zone_ ) {};
  void operator()( std::ostream& out ) const
  {
    out.precision(4);
    out << "graph [bgcolor=lightgrey];" << std::endl;
    out <<
      boost::format( "node [shape=box color=\"%s\"];" ) % get_bounding_color()
      << std::endl;
    out << "label = \"latex\";" << std::endl;
    out <<
      boost::format(
        "texlbl = \"\\huge{$time(s) = %g \
         \\ \\ \\ \\ T_9 = %g \
         \\ \\ \\ \\ \\rho(g/cc) = %g \
         \\ \\ \\ \\ {\\mathrm{flow}_{max}} = %g$}\";"
      ) %
      boost::lexical_cast<double>( zone.getProperty( nnt::s_TIME ) ) %
      boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) %
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) %
      boost::lexical_cast<double>( zone.getProperty( "max flow" ) ) <<
      std::endl;
  };
  nnt::Zone& zone;
};

//##############################################################################
// my_output().
//##############################################################################

void
my_output(
  my_graph_t &g,
  nnt::Zone zone,
  const char * s_output
)
{

  std::ofstream my_out;
  my_out.open( s_output );
  boost::write_graphviz(
    my_out,
    g,
    my_vertex_writer( g ),
    my_edge_writer( g ),
    my_graph_writer( zone )
  );
  my_out.close();

}

//##############################################################################
// create_graph().
//##############################################################################

my_graph_t
create_graph(
  nnt::Zone& zone,
  Libnucnet__NetView * p_view,
  std::string s_type
)
{

  my_graph_t g;

  my_graph_t::vertex_descriptor v;

  my_vertex_hash_t vertices_hash;

  nnt::species_list_t species_list =
    nnt::make_species_list(
      Libnucnet__Net__getNuc(
        Libnucnet__NetView__getNet( p_view )
      )
    );

  BOOST_FOREACH( nnt::Species species, species_list )
  {

    v = boost::add_vertex( g );

    g[v].setNucnetSpecies( species.getNucnetSpecies() );

    vertices_hash.insert(
      my_graph_vertex_t(
        Libnucnet__Species__getName( species.getNucnetSpecies() ),
        species.getNucnetSpecies(),
        v
      )
    );

  }

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list(
      Libnucnet__Net__getReac(
        Libnucnet__NetView__getNet( p_view )
      )
    );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    std::pair<double,double> flows =
      user::compute_flows_for_reaction(
        zone,
        reaction.getNucnetReaction()
      );

    if( s_type == "net" )
      add_reaction_net_flow_edges( g, vertices_hash, reaction, flows );
    else if( s_type == "all" )
      add_reaction_all_flow_edges( g, vertices_hash, reaction, flows );
    else
    {
      std::cerr << "No such flow type." << std::endl;
      exit( EXIT_FAILURE );
    }

  }

  return g;

}

//##############################################################################
// main().
//##############################################################################

int
main( int argc, char **argv )
{

  Libnucnet *p_my_nucnet;
  Libnucnet__NetView * p_view;
  Libnucnet__NucView * p_nuc_view;

  my_graph_t g, sub_g;

  if( argc != 9 )
  {
    fprintf(
      stderr,
      "\nUsage: %s xml_file nuc_xpath reac_xpath zone_xpath induced_nuc_xpath scaling flow_type dot_file_base\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  xml_file = network xml file\n\n"
    );
    fprintf(
      stderr,
      "  nuc_xpath = XPath expression to select nuclides\n\n"
    );
    fprintf(
      stderr,
      "  reac_xpath = XPath expression to select reactions\n\n"
    );
    fprintf(
      stderr,
      "  zone_xpath = XPath expression to select zones\n\n"
    );
    fprintf(
      stderr,
      "  induced_nuc_xpath = XPath expression to induced subgraph\n\n"
    );
    fprintf(
      stderr,
      "  flow_type = flow type (\"net\" or \"all\")\n\n"
    );
    fprintf(
      stderr,
      "  scaling = scaling of weights (\"linear\" or \"logarithmic\")\n\n"
    );
    fprintf(
      stderr,
      "  dot_file_base = base name for output dot file\n\n"
    );
    return EXIT_FAILURE;
  }

  p_my_nucnet =
    Libnucnet__new_from_xml(
      argv[1],
      NULL,
      NULL,
      argv[4]
    );

  //============================================================================
  // Get net view.
  //============================================================================
  
  p_view =
    Libnucnet__NetView__new(
      Libnucnet__getNet( p_my_nucnet ),
      argv[2],
      argv[3]
    );

  //============================================================================
  // Create reaction color map.
  //============================================================================
  
  create_reaction_color_map(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

  //============================================================================
  // Create reaction linestyle map.
  //============================================================================

  create_reaction_linestyle_map(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

  //============================================================================
  // Register rate functions.  Set two-d weak rates hash.
  //============================================================================
  
  user::register_my_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

#ifdef MY_USER
  my_user::register_my_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );
#endif
  
  user::set_two_d_weak_rates_hashes(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

  //============================================================================
  // Get subgraph view.
  //============================================================================
  
  p_nuc_view =
    Libnucnet__NucView__new(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
      argv[5]
    );

  //============================================================================
  // Get zone list.
  //============================================================================
  
  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  nnt::zone_list_t zone_list = nnt::make_zone_list( p_my_nucnet );

  //============================================================================
  // Iterate zones.
  //============================================================================
  
  nnt::Zone zone_first = *(zone_list.begin());

  BOOST_FOREACH( nnt::Zone zone, zone_list )
  {

    std::cout << Libnucnet__Zone__getLabel( zone.getNucnetZone(), 1 ) << std::endl;

    if( !( zone == zone_first ) )
    {
      Libnucnet__Zone__copy_net_views(
        zone.getNucnetZone(),
        zone_first.getNucnetZone()
      );
    }

    user::update_my_rate_functions_data( zone );

    g = create_graph( zone, p_view, argv[6] );

    sub_g = create_induced_subgraph( g, p_nuc_view, true );

    //==========================================================================
    // Record maximum flow and scale subgraph weights.
    //==========================================================================

    zone.updateProperty(
      "max flow",
      boost::lexical_cast<std::string>(
        get_max_flow( sub_g )
      )
    );

    scale_graph_weights( sub_g, argv[7] );

    //==========================================================================
    // Induce another subgraph to get rid of isolated vertices, if desired.
    // Add anchors to fix graph boundaries.
    //==========================================================================
    
    if( !B_ALLOW_ISOLATED_VERTICES )
    {
      sub_g =
        create_induced_subgraph( sub_g, p_nuc_view, false );
    }

    //==========================================================================
    // Add solar coloring.
    //==========================================================================

    add_solar( sub_g, p_nuc_view );
    
    //==========================================================================
    // Output.
    //==========================================================================
    
    std::stringstream ss;

    ss <<
       std::setw(5) <<
       std::setfill('0') <<
       Libnucnet__Zone__getLabel( zone.getNucnetZone(), 1 );

    std::string s_label(ss.str());

    std::string s_output =
      std::string( argv[8] ) +
      std::string( "_" ) +
      s_label;
   
    my_output( sub_g, zone, s_output.c_str() );

  }

  //============================================================================
  // Clean up.
  //============================================================================

  Libnucnet__NucView__free( p_nuc_view );
  Libnucnet__NetView__free( p_view );
  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}

