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
//! \brief Example code to write out the integrated current net flow graph
//!    of a calculation in dot format.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Include.
//##############################################################################

#include <fstream>
#include <ostream>
#include <iostream>
#include <sstream>
#include <iomanip>

#include <map>

#include <boost/graph/graphviz.hpp>

#include "nnt/iter.h"

#include "graph_types.h"
#include "graph_helper.h"

#include "user/user_rate_functions.h"

#define D_SCALE  100
#define D_MAX_FLOW_SCALE 10

#define B_ALLOW_ISOLATED_VERTICES false

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
// fill_flow_current_hash().
//##############################################################################

void
fill_flow_current_hash(
  const char * s_name,
  const char * s_tag1,
  const char * s_tag2,
  const char * s_value,
  current_multi& currents
)
{

  currents.insert(
    Reaction_Current(
      s_tag1,
      boost::lexical_cast<double>( s_value )
    )   
  );

}

//##############################################################################
// create_graph().
//##############################################################################

my_graph_t
create_graph(
  Libnucnet__NetView * p_view,
  current_multi currents
)
{

  my_graph_t g;

  my_graph_t::vertex_descriptor v;
  my_graph_t::edge_descriptor e;

  bool b_added_edge;

  my_vertex_hash_t vertices_hash;
  boost::unordered_set<Graph_Link> my_links;

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

  BOOST_FOREACH( Reaction_Current my_current, currents )
  {

    Libnucnet__Reaction * p_reaction =
      Libnucnet__Reac__getReactionByString(
        Libnucnet__Net__getReac(
          Libnucnet__NetView__getNet( p_view )
        ),
        my_current.name.c_str()
      );

    if( p_reaction )
    {

      if( my_current.current >= 0 )
      {

	nnt::reaction_element_list_t reactant_list =
	  nnt::make_reaction_nuclide_reactant_list( p_reaction );

	BOOST_FOREACH( nnt::ReactionElement reactant, reactant_list )
	{

	  my_vertex_hash_t::iterator itr =
	    vertices_hash.find(
	      Libnucnet__Reaction__Element__getName(
		reactant.getNucnetReactionElement()
	      )
	    );

	  nnt::reaction_element_list_t product_list =
	    nnt::make_reaction_nuclide_product_list( p_reaction );

	  BOOST_FOREACH( nnt::ReactionElement product, product_list )
	  {

	    my_vertex_hash_t::iterator itp =
	      vertices_hash.find(
		Libnucnet__Reaction__Element__getName(
		  product.getNucnetReactionElement()
		)
	      );

	    boost::unordered_set<Graph_Link>::iterator itl =
	      my_links.find(
		Graph_Link(
		  Libnucnet__Reaction__Element__getName(
		    reactant.getNucnetReactionElement()
		  ),
		  Libnucnet__Reaction__Element__getName(
		    product.getNucnetReactionElement()
		  )
		)
	      );

	    if( itl == my_links.end() )
	    {

	      boost::tie( e, b_added_edge ) =
		boost::add_edge( itr->v, itp->v, g );
	      g[e].setNucnetReaction( p_reaction );
	      g[e].setWeight( my_current.current );

	      my_links.insert(
		Graph_Link(
		  Libnucnet__Reaction__Element__getName(
		    reactant.getNucnetReactionElement()
		  ),
		  Libnucnet__Reaction__Element__getName(
		    product.getNucnetReactionElement()
		  )
		)
	      );

	    }

	  }

	}

      }
      else
      {

	nnt::reaction_element_list_t product_list =
	  nnt::make_reaction_nuclide_product_list( p_reaction );

	BOOST_FOREACH( nnt::ReactionElement product, product_list )
	{

	  my_vertex_hash_t::iterator itp =
	    vertices_hash.find(
	      Libnucnet__Reaction__Element__getName(
		product.getNucnetReactionElement()
	      )
	    );

	  nnt::reaction_element_list_t reactant_list =
	    nnt::make_reaction_nuclide_reactant_list( p_reaction );

	  BOOST_FOREACH( nnt::ReactionElement reactant, reactant_list )
	  {

	    my_vertex_hash_t::iterator itr =
	      vertices_hash.find(
		Libnucnet__Reaction__Element__getName(
		  reactant.getNucnetReactionElement()
		)
	      );

	    boost::unordered_set<Graph_Link>::iterator itl =
	      my_links.find(
		Graph_Link(
		  Libnucnet__Reaction__Element__getName(
		    product.getNucnetReactionElement()
		  ),
		  Libnucnet__Reaction__Element__getName(
		    reactant.getNucnetReactionElement()
		  )
		)
	      );

	    if( itl == my_links.end() )
	    {

	      boost::tie( e, b_added_edge ) =
		boost::add_edge( itp->v, itr->v, g );
	      g[e].setNucnetReaction( p_reaction );
	      g[e].setWeight( -my_current.current );

	      my_links.insert(
		Graph_Link(
		  Libnucnet__Reaction__Element__getName(
		    product.getNucnetReactionElement()
		  ),
		  Libnucnet__Reaction__Element__getName(
		    reactant.getNucnetReactionElement()
		  )
		)
	      );

	    }

	  }

	}

      }

    }

    my_links.clear();

  }

  return g;

}

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
  void operator()( std::ostream& out ) const
  {
    out << "graph [bgcolor=lightgrey];" << std::endl;
    out << "node [shape=box color=" << get_bounding_color() <<
           "];" << std::endl;
    out.precision(4);
  };
};

//##############################################################################
// my_output().
//##############################################################################

void
my_output(
  my_graph_t &g,
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
    my_graph_writer( )
  );
  my_out.close();

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
  Libnucnet__Zone * p_flow_current_zone;
  current_multi currents;

  my_graph_t g;

  if( argc != 7 )
  {
    fprintf(
      stderr,
      "\nUsage: %s xml_file nuc_xpath reac_xpath induced_nuc_xpath scaling dot_file_base\n\n",
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
      "  induced_nuc_xpath = XPath expression to induced subgraph\n\n"
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
      NULL
    );

  p_view =
    Libnucnet__NetView__new(
      Libnucnet__getNet( p_my_nucnet ),
      argv[2],
      argv[3]
    );

  //============================================================================
  // Get currents zone.
  //============================================================================
  
  p_flow_current_zone =
    Libnucnet__getZoneByLabels( p_my_nucnet, "integrated currents", "0", "0" );

  //============================================================================
  // Fill currents hash.
  //============================================================================

  Libnucnet__Zone__iterateOptionalProperties(
    p_flow_current_zone,
    nnt::s_FLOW_CURRENT,
    NULL,
    NULL,
    (Libnucnet__Zone__optional_property_iterate_function)
      fill_flow_current_hash,
    &currents
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
  // Printout sorted currents.
  //============================================================================

  BOOST_FOREACH( Reaction_Current my_current, currents )
  {

    if(
      Libnucnet__Reac__getReactionByString(
        Libnucnet__Net__getReac( Libnucnet__NetView__getNet( p_view ) ),
        my_current.name.c_str()
      )
    )
    {

      if( fabs( my_current.current ) > 1.e-20 )
      {
        std::cout <<
          boost::format( "%-60s  %.5e\n" ) %
            my_current.name %
            my_current.current;
      }

    }

  }
      
  //============================================================================
  // Get subgraph view.
  //============================================================================
  
  g = create_graph( p_view, currents );

  p_nuc_view =
    Libnucnet__NucView__new(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
      argv[4]
    );

  my_graph_t sub_g = create_induced_subgraph( g, p_nuc_view, true );

  //============================================================================
  // Record maximum flow and scale subgraph weights.
  //============================================================================

  scale_graph_weights( sub_g, argv[5] );

  //============================================================================
  // Induce another subgraph to get rid of isolated vertices, if desired.
  // Add anchors to fix graph boundaries.
  //============================================================================
    
  if( !B_ALLOW_ISOLATED_VERTICES )
  {
    sub_g = create_induced_subgraph( sub_g, p_nuc_view, false );
  }

  //============================================================================
  // Add solar coloring.
  //============================================================================

  add_solar( sub_g, p_nuc_view );
    
  //============================================================================
  // Output.
  //============================================================================
    
  my_output( sub_g, argv[6] );

  //============================================================================
  // Clean up.
  //============================================================================

  Libnucnet__NucView__free( p_nuc_view );
  Libnucnet__NetView__free( p_view );
  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
