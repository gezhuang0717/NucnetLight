///////////////////////////////////////////////////////////////////////////////
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
//! \brief Example code to write out the abundance graph of a zone in dot
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

#include "graph_helper.h"

#define D_BASE  120.
#define D_SCALE 100

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
         \\ \\ \\ \\ \\rho(g/cc) = %g$}\";"
      ) %
      boost::lexical_cast<double>( zone.getProperty( nnt::s_TIME ) ) %
      boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) %
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) <<
      std::endl;
  };
  nnt::Zone& zone;
};

//##############################################################################
// my_output().
//##############################################################################

void
my_output(
  my_graph_t& g,
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
    boost::default_writer(),
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
  Libnucnet__NucView * p_nuc_view,
  double d_cutoff
)
{

  std::stringstream s_base;
  std::string s_color, s_white = "#FFFFFF";
  my_graph_t g;

  my_graph_t::vertex_descriptor v;

  nnt::species_list_t species_list =
    nnt::make_species_list(
      Libnucnet__NucView__getNuc( p_nuc_view )
    );

  BOOST_FOREACH( nnt::Species species, species_list )
  {

    v = boost::add_vertex( g );

    g[v].setNucnetSpecies( species.getNucnetSpecies() );

    double d_y =
      Libnucnet__Zone__getSpeciesAbundance(
        zone.getNucnetZone(),
        species.getNucnetSpecies()
      );

    if( d_y < d_cutoff )
    {
      s_color = s_white;
    }
    else
    {
      s_base.str("");
      s_base << std::hex <<
        static_cast<int>(
          D_BASE + (256. - D_BASE) *
          log10( d_y ) / log10( d_cutoff )
        );
      s_color = "#" + s_base.str() + s_base.str() + s_base.str();
    }

    g[v].updateExtraData( my_vertex_data( s_color, "box" ) );

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
  Libnucnet__NucView * p_nuc_view;

  my_graph_t g;

  if( argc != 6 )
  {
    fprintf(
      stderr,
      "\nUsage: %s xml_file nuc_xpath zone_xpath cutoff dot_file_base\n\n",
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
      "  zone_xpath = XPath expression to select zones\n\n"
    );
    fprintf(
      stderr,
      "  cutoff = minimum abundance to include\n\n"
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
      argv[3]
    );

  //============================================================================
  // Get subgraph view.
  //============================================================================
  
  p_nuc_view =
    Libnucnet__NucView__new(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
      argv[2]
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
  
  BOOST_FOREACH( nnt::Zone zone, zone_list )
  {

    g = create_graph( zone, p_nuc_view, atof( argv[4] ) );

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
      std::string( argv[5] ) +
      std::string( "_" ) +
      s_label;
   
    my_output( g, zone, s_output.c_str() );

  }

  //============================================================================
  // Clean up.
  //============================================================================

  Libnucnet__NucView__free( p_nuc_view );
  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
