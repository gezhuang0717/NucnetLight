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
//! \brief Header to create graphs from zone files.
////////////////////////////////////////////////////////////////////////////////

#ifndef GRAPH_H
#define GRAPH_H

//##############################################################################
// Includes.
//##############################################################################

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/scoped_ptr.hpp>

#include "nnt/auxiliary.h"
#include "nnt/iter.h"

#include <Libnucnet.h>

namespace nnt
{

template<typename T>
class Graph_Species
{
  public:
    Graph_Species(){};
    ~Graph_Species(){};
    void setNucnetSpecies( Libnucnet__Species * p ){ pSpecies = p; };
    void updateExtraData( T x )
    {
      p_x.reset( new T );
      *p_x = x;
    }
    Libnucnet__Species * getNucnetSpecies() const { return pSpecies; }
    T getExtraData() const { return *p_x; }

  private:
    Libnucnet__Species * pSpecies;
    boost::shared_ptr<T> p_x;

};

template<typename T>
class Graph_Reaction {
  public:
    Graph_Reaction(){};
    ~Graph_Reaction(){};
    void setNucnetReaction( Libnucnet__Reaction * p ){ pReaction = p; };
    void setWeight( double x ){ dWeight = x; } ;
    void setExtraData( T x )
    {
      p_x.reset( new T );
      *p_x = x;
    }
    Libnucnet__Reaction * getNucnetReaction() const { return pReaction; }
    double getWeight() const { return dWeight; }
    T getExtraData() const { return *p_x; }

  private:
    Libnucnet__Reaction * pReaction;
    double dWeight;
    boost::shared_ptr<T> p_x;

};

/**
 * A Boost graph type representing a network graph.
 */

template<typename T1, typename T2>
class net_graph_base
{
  public:
    typedef boost::adjacency_list<
      boost::listS,
      boost::vecS,
      boost::bidirectionalS,
      Graph_Species<T1>,
      Graph_Reaction<T2>
    > net_graph;
};

template<typename T1, typename T2>
class net_graph
{
  public:
    typedef typename net_graph_base<T1,T2>::net_graph type;
};

template<typename T1, typename T2>
struct Graph_Vertex_t
{
  std::string name;
  Libnucnet__Species * pSpecies;
  typename net_graph<T1,T2>::type::vertex_descriptor v;

  Graph_Vertex_t(
    const std::string &s,
    Libnucnet__Species * p_s,
    typename net_graph<T1,T2>::type::vertex_descriptor v_in
  ) : name( s ), pSpecies( p_s ), v( v_in ){}
};

template<typename T1, typename T2>
struct Graph_Vertex
{
  typedef Graph_Vertex_t<T1,T2> type;
};
  
template<typename T1, typename T2>
struct vertex_multi_base
{
  typedef
    boost::multi_index::multi_index_container<
      typename Graph_Vertex<T1,T2>::type,
      boost::multi_index::indexed_by<
	boost::multi_index::hashed_unique<
	  boost::multi_index::member< 
	    typename Graph_Vertex<T1,T2>::type,
	    std::string,
	    &Graph_Vertex<T1,T2>::type::name
	  > 
	>,
	boost::multi_index::hashed_unique< 
	  boost::multi_index::member<  
	    typename Graph_Vertex<T1,T2>::type,
	    Libnucnet__Species *,
	    &Graph_Vertex<T1,T2>::type::pSpecies
	  > 
	>, 
	boost::multi_index::hashed_unique< 
	  boost::multi_index::member<  
	    typename Graph_Vertex<T1,T2>::type,
	    typename net_graph<T1,T2>::type::vertex_descriptor,
	    &Graph_Vertex<T1,T2>::type::v
	  > 
	> 
      > 
    > vertex_multi;
};

template<typename T1, typename T2>
class vertex_multi
{
  public:
    typedef typename vertex_multi_base<T1,T2>::vertex_multi type;
};

} // namespace nnt

#endif // GRAPH_H
