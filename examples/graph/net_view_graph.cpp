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
//! \brief Example code to write out a graph of a network view in dot
//!    format.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Include.
//##############################################################################

#include <fstream>
#include <ostream>
#include <iostream>

#include <boost/graph/graphviz.hpp>
#include "graph_helper.h"

#define D_SCALE  75 

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
  my_edge_writer(
    my_graph_t& g_,
    Libnucnet__Net * p_net_
  ) : g( g_ ), pNet( p_net_ ) {};
  template <class Edge>
  void operator()( std::ostream& out, Edge e )
  {
    std::string s_color = get_reaction_color( g[e].getNucnetReaction() );
    char * s_latex_string =
      Libnucnet__Net__createValidReactionLatexString(
        pNet,
        g[e].getNucnetReaction()
      );
    out <<
      boost::format(
        " [ style=\"line width = 1.0pt\", color = \"%s\" ]"
      ) %
        s_color
      << std::endl;
    free( s_latex_string );
  }
  my_graph_t& g;
  Libnucnet__Net *pNet;
};

//##############################################################################
// The graph writer.
//##############################################################################

struct
my_graph_writer
{
  void operator()( std::ostream& out ) const
  {
    out << "graph [bgcolor=lightgrey]" << std::endl;
    out << "node [shape=box color=blue]" << std::endl;
  };
} my_graph_writer;

//##############################################################################
// my_output().
//##############################################################################

void
my_output(
  my_graph_t& g,
  Libnucnet__Net * p_net,
  const char * s_output
)
{

  std::ofstream my_out;
  my_out.open( s_output );
  boost::write_graphviz(
    my_out,
    g,
    my_vertex_writer( g ),
    my_edge_writer( g, p_net ),
/*
    boost::default_writer(),
*/
    my_graph_writer
  );
  my_out.close();

}

//##############################################################################
// create_net_view_graph().
//##############################################################################

my_graph_t
create_net_view_graph( Libnucnet__NetView * p_view )
{

  my_vertex_hash_t vertices_hash;

  my_graph_t g;

  my_graph_t::vertex_descriptor v;
  my_graph_t::edge_descriptor e;

  bool b_added_edge;

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

    nnt::reaction_element_list_t reactant_list =
      nnt::make_reaction_nuclide_reactant_list(
        reaction.getNucnetReaction()
      );

    BOOST_FOREACH( nnt::ReactionElement reactant, reactant_list )
    {

      my_vertex_hash_t::iterator itr =
	vertices_hash.find(
	  Libnucnet__Reaction__Element__getName(
	    reactant.getNucnetReactionElement()
	  )
	);

      nnt::reaction_element_list_t product_list =
        nnt::make_reaction_nuclide_product_list(
          reaction.getNucnetReaction()
        );

      BOOST_FOREACH( nnt::ReactionElement product, product_list )
      {

	my_vertex_hash_t::iterator itp =
	  vertices_hash.find(
	    Libnucnet__Reaction__Element__getName(
	      product.getNucnetReactionElement()
	    )
	  );

	boost::tie( e, b_added_edge ) = boost::add_edge( itr->v, itp->v, g );
        g[e].setNucnetReaction( reaction.getNucnetReaction() );

      }

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

  Libnucnet__Net *p_my_net;
  Libnucnet__NetView * p_view;
  Libnucnet__NucView * p_nuc_view;

  my_graph_t g;

  my_graph_t::edge_iterator
    ei, ei_end, e_next;

  if( argc < 5 || argc > 6 )
  {
    fprintf(
      stderr,
      "\nUsage: %s xml_file nuc_xpath reac_xpath dot_file\n\n",
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
      "  induced_nuc_xpath = XPath expression to induced subgraph (optional)\n\n"
    );
    fprintf(
      stderr,
      "  dot_file = output dot file\n\n"
    );
    return EXIT_FAILURE;
  }

  p_my_net =
    Libnucnet__Net__new_from_xml(
      argv[1],
      NULL,
      NULL
    );

  p_view =
    Libnucnet__NetView__new(
      p_my_net,
      argv[2],
      argv[3]
    );

  //============================================================================
  // Create reaction color map.
  //============================================================================
  
  create_reaction_color_map(
    Libnucnet__Net__getReac( Libnucnet__NetView__getNet( p_view ) )
  );

  g = create_net_view_graph( p_view );

  if( argc == 5 )
    my_output( g, p_my_net, argv[4] );
  else
  {
    p_nuc_view =
      Libnucnet__NucView__new(
        Libnucnet__Net__getNuc( Libnucnet__NetView__getNet( p_view ) ), 
        argv[4]
      );
    my_graph_t sub_g =
      create_induced_subgraph( g, p_nuc_view, true );
    add_solar( sub_g, p_nuc_view );
    Libnucnet__NucView__free( p_nuc_view );
    my_output( sub_g, p_my_net, argv[5] );
  }
 
  //============================================================================
  // Clean up.
  //============================================================================

  Libnucnet__NetView__free( p_view );
  Libnucnet__Net__free( p_my_net );

  return EXIT_SUCCESS;

}
