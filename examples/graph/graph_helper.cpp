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
//! \brief Helper file for graph routines.
////////////////////////////////////////////////////////////////////////////////

#include "graph_helper.h"

my_vertex_data
set_vertex_data( std::string s_color )
{
  return my_vertex_data( s_color, "box" );

} 

//##############################################################################
// add_solar().
//##############################################################################

void
add_solar(
  my_graph_t &g,
  Libnucnet__NucView * p_nuc_view
)
{
  
  my_graph_t::vertex_iterator vi, vi_end;
  my_graph_t::vertex_descriptor v;
  Libnucnet__Species * p_species;
  boost::unordered_map<std::string,my_graph_t::vertex_descriptor>
    my_vertices;
  boost::unordered_map<std::string,std::string> special_vertex_color_map;
  boost::unordered_map<std::string,std::string>::iterator its;
  std::string s_shape = "box";

  std::pair<std::string,std::string> colors = solar_color();

  for( boost::tie( vi, vi_end ) = vertices( g ); vi != vi_end; vi++ )
  {
    my_vertices[Libnucnet__Species__getName( g[*vi].getNucnetSpecies() )] = *vi;
 //   g[*vi].updateExtraData( set_vertex_data( colors.first ) );
      g[*vi].updateExtraData( my_vertex_data( colors.first, "box" ) );
  }

  std::vector<std::string> stables = nnt::get_stable_species();

  BOOST_FOREACH( std::string s, stables )
  {

    if( my_vertices.find( s ) != my_vertices.end() )
    {
      g[my_vertices[ s ]].updateExtraData( set_vertex_data( colors.second ) );
    }
    else
    {

      p_species =
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__NucView__getNuc( p_nuc_view ),
          s.c_str()
        );
    
      if( p_species )
      {

        v = add_vertex( g );
        g[v].setNucnetSpecies( p_species );
        g[v].updateExtraData( set_vertex_data( colors.second ) );
        my_vertices[Libnucnet__Species__getName( p_species )] = v;

      }

    }
        
  }

  special_vertex_color_map = get_special_vertex_color_map();

  for(
    its = special_vertex_color_map.begin();
    its != special_vertex_color_map.end();
    its++
  )
  {

    if( my_vertices.find( its->first ) != my_vertices.end() )
    {
      g[my_vertices[its->first]].updateExtraData(
        set_vertex_data( its->second )
      );
    }

  }

}

//##############################################################################
// hash_value().
//##############################################################################

std::size_t hash_value(Graph_Link const &p)
{
  std::size_t seed = 0;
  boost::hash_combine( seed, p.out_name );
  boost::hash_combine( seed, p.in_name );
  return seed;
}

//##############################################################################
// create_induced_subgraph().
//##############################################################################

/**
 * \brief Returns an induced subgraph of the input graph.
 * \param g The input graph.
 * \param p_view A Libnucnet__NucView giving the vertices of the induced
 *                 subgraph.
 * \param b_allow_isolated A boolean to allow isolated vertices in the subgraph
 *                         (true) or not (false).
 * \return A net_graph that is the induced subgraph of the input graph.
 */

my_graph_t
create_induced_subgraph(
  my_graph_t& g,
  Libnucnet__NucView * p_view,
  bool b_allow_isolated
)
{

  my_graph_t sub_g;

  my_graph_t::vertex_descriptor v, v1, v2;
  my_graph_t::vertex_iterator vi, vi_end;

  my_graph_t::edge_descriptor e;
  my_graph_t::edge_iterator ei, ei_end;

  bool b_added_edge;

  my_vertex_hash_t vertices_hash;

  for( boost::tie( vi, vi_end ) = boost::vertices( g ); vi != vi_end; vi++ )
  {

    if(
      Libnucnet__Nuc__getSpeciesByName(
        Libnucnet__NucView__getNuc( p_view ),
        Libnucnet__Species__getName( g[*vi].getNucnetSpecies() )
      )
    )
    {

      if( 
        b_allow_isolated ||
        (
          boost::in_degree( *vi, g ) != 0 ||
          boost::out_degree( *vi, g ) != 0
        )
      )
      {

        v = boost::add_vertex( sub_g );

        sub_g[v].setNucnetSpecies( g[*vi].getNucnetSpecies() );

        vertices_hash.insert(
          my_graph_vertex_t(
            Libnucnet__Species__getName( g[*vi].getNucnetSpecies() ),
            g[*vi].getNucnetSpecies(),
            v
          )
        );

      }

    }

  }

  for( tie( ei, ei_end ) = boost::edges( g ); ei != ei_end; ei++ )
  {

    v1 = boost::source( *ei, g );
    v2 = boost::target( *ei, g );

    if(
      Libnucnet__Nuc__getSpeciesByName(
        Libnucnet__NucView__getNuc( p_view ),
        Libnucnet__Species__getName( g[v1].getNucnetSpecies() )
      ) &&
      Libnucnet__Nuc__getSpeciesByName(
        Libnucnet__NucView__getNuc( p_view ),
        Libnucnet__Species__getName( g[v2].getNucnetSpecies() )
      )
    )
    {

      my_vertex_hash_t::iterator its =
	vertices_hash.find(
          Libnucnet__Species__getName( g[v1].getNucnetSpecies() )
	);

      my_vertex_hash_t::iterator itt =
	vertices_hash.find(
          Libnucnet__Species__getName( g[v2].getNucnetSpecies() )
	);

      tie( e, b_added_edge ) = boost::add_edge( its->v, itt->v, sub_g );

      g[e].setNucnetReaction( g[*ei].getNucnetReaction() );
      g[e].setWeight( g[*ei].getWeight() );

    }

  }

  //============================================================================
  // Add anchors if not already present.
  //============================================================================
  
  std::vector<Libnucnet__Species *>
  my_vector =
    nnt::get_nuc_anchors(
      Libnucnet__NucView__getNuc( p_view )
    );

  BOOST_FOREACH( Libnucnet__Species * p_species, my_vector )
  {

    my_vertex_hash_t::iterator it =
      vertices_hash.find(
        Libnucnet__Species__getName( p_species )
      );

    if( it == vertices_hash.end() )
    {
      v = add_vertex( sub_g );
      sub_g[v].setNucnetSpecies( p_species );
      vertices_hash.insert(
        my_graph_vertex_t(
          Libnucnet__Species__getName( p_species ),
          g[v].getNucnetSpecies(),
          v
        )
      );
    }

  }

  //============================================================================
  // Return.
  //============================================================================
  
  return sub_g;

}

//##############################################################################
// get_max_flow().
//##############################################################################

double
get_max_flow(
  my_graph_t& g
)
{

  my_graph_t::edge_iterator ei, ei_end;

  double d_max = 0.;

  for(
    boost::tie( ei, ei_end ) = boost::edges( g );
    ei != ei_end;
    ei++
  )
  {
    if( fabs( g[*ei].getWeight() ) > d_max ) d_max =
      fabs( g[*ei].getWeight() );
  }

  return d_max;

}

//##############################################################################
// add_reaction_net_flow_edges().
//##############################################################################

void
add_reaction_net_flow_edges(
  my_graph_t& g,
  my_vertex_hash_t& vertices_hash,
  nnt::Reaction& reaction,
  std::pair<double,double>& flows
)
{

  my_graph_t::edge_descriptor e;

  bool b_added_edge;
  boost::unordered_set<Graph_Link> my_links;

  double d_net_flow = flows.first - flows.second;

  if( d_net_flow >= 0 )
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
	  g[e].setNucnetReaction( reaction.getNucnetReaction() );
	  g[e].setWeight( d_net_flow );

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
	  g[e].setNucnetReaction( reaction.getNucnetReaction() );
	  g[e].setWeight( -d_net_flow );

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

//##############################################################################
// add_reaction_all_flow_edges().
//##############################################################################

void
add_reaction_all_flow_edges(
  my_graph_t& g,
  my_vertex_hash_t& vertices_hash,
  nnt::Reaction& reaction,
  std::pair<double,double>& flows
)
{

  my_graph_t::edge_descriptor e;

  bool b_added_edge;
  boost::unordered_set<Graph_Link> my_links;
  boost::unordered_set<Graph_Link>::iterator itl;

  my_vertex_hash_t::iterator itr;
  my_vertex_hash_t::iterator itp;

  nnt::reaction_element_list_t reactant_list =
    nnt::make_reaction_nuclide_reactant_list(
      reaction.getNucnetReaction()
    );

  nnt::reaction_element_list_t product_list =
    nnt::make_reaction_nuclide_product_list(
      reaction.getNucnetReaction()
    );

  BOOST_FOREACH( nnt::ReactionElement reactant, reactant_list )
  {

    itr =
      vertices_hash.find(
	Libnucnet__Reaction__Element__getName(
	  reactant.getNucnetReactionElement()
	)
      );

    BOOST_FOREACH( nnt::ReactionElement product, product_list )
    {

      itp =
        vertices_hash.find(
          Libnucnet__Reaction__Element__getName(
            product.getNucnetReactionElement()
          )
        );

      itl =
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

	boost::tie( e, b_added_edge ) = boost::add_edge( itr->v, itp->v, g );
	g[e].setNucnetReaction( reaction.getNucnetReaction() );
	g[e].setWeight( flows.first );

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

  BOOST_FOREACH( nnt::ReactionElement product, product_list )
  {

    itp =
      vertices_hash.find(
	Libnucnet__Reaction__Element__getName(
	  product.getNucnetReactionElement()
	)
      );

    BOOST_FOREACH( nnt::ReactionElement reactant, reactant_list )
    {

      itr =
        vertices_hash.find(
          Libnucnet__Reaction__Element__getName(
            reactant.getNucnetReactionElement()
          )
        );

      itl =
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

	boost::tie( e, b_added_edge ) = boost::add_edge( itp->v, itr->v, g );
	g[e].setNucnetReaction( reaction.getNucnetReaction() );
	g[e].setWeight( flows.second );

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

