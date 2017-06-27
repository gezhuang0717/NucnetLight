///////////////////////////////////////////////////////////////////////////////
// This file was originally written by Bradley S. Meyer.
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
////////////////////////////////////////////////////////////////////////////////
 
////////////////////////////////////////////////////////////////////////////////
//! \file
//! \brief Example code to mix mass and species over zones.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <boost/bind.hpp>

#include "nnt/write_output_xml.h"

#include "user/hdf5_routines.h"
#include "user/multi_zone_utilities.h"

//##############################################################################
// Defines.
//##############################################################################

#define D_MIX_RATE  1.e-02

//##############################################################################
// Types.
//##############################################################################

typedef
boost::tuple< std::string, std::string, std::string, double > link_tuple_t;

typedef std::map<size_t, link_tuple_t> link_map_t;

//##############################################################################
// set_iterative_solver_params().
//##############################################################################

boost::unordered_map<std::string,std::string>
set_iterative_solver_parameters( std::vector<nnt::Zone>& zones )
{

  boost::unordered_map<std::string,std::string> param_map;

  //============================================================================
  // Set parameters for iterative solver.
  //============================================================================

  if( zones[0].hasProperty( nnt::s_ITER_SOLVER ) )
  {
    param_map.insert(
      boost::unordered_map<std::string,std::string>::value_type(
	nnt::s_ITER_SOLVER,
	zones[0].getProperty( nnt::s_ITER_SOLVER )
      )
    );
  }
  else
  {
    param_map.insert(
      boost::unordered_map<std::string,std::string>::value_type(
	nnt::s_ITER_SOLVER,
	"gmres"
      )
    );
  }

  if( zones[0].hasProperty( nnt::s_ITER_SOLVER_REL_TOL ) )
  {
    param_map.insert(
      boost::unordered_map<std::string,std::string>::value_type(
	nnt::s_ITER_SOLVER_REL_TOL,
	zones[0].getProperty( nnt::s_ITER_SOLVER_REL_TOL )
      )
    );
  }
  else
  {
    param_map.insert(
      boost::unordered_map<std::string,std::string>::value_type(
	nnt::s_ITER_SOLVER_REL_TOL,
	"1.e-12"
      )
    );
  }

  if( zones[0].hasProperty( nnt::s_ITER_SOLVER_ABS_TOL ) )
  {
    param_map.insert(
      boost::unordered_map<std::string,std::string>::value_type(
	nnt::s_ITER_SOLVER_ABS_TOL,
	zones[0].getProperty( nnt::s_ITER_SOLVER_ABS_TOL )
      )
    );
  }
  else
  {
    param_map.insert(
      boost::unordered_map<std::string,std::string>::value_type(
	nnt::s_ITER_SOLVER_ABS_TOL,
	"1.e-25"
      )
    );
  }

  if( zones[0].hasProperty( nnt::s_ITER_SOLVER_MAX_ITERATIONS ) )
  {
    param_map.insert(
      boost::unordered_map<std::string,std::string>::value_type(
	nnt::s_ITER_SOLVER_MAX_ITERATIONS,
	zones[0].getProperty( nnt::s_ITER_SOLVER_MAX_ITERATIONS )
      )
    );
  }
  else
  {
    param_map.insert(
      boost::unordered_map<std::string,std::string>::value_type(
	nnt::s_ITER_SOLVER_MAX_ITERATIONS,
	"50"
      )
    );
  }

  if( zones[0].hasProperty( nnt::s_ILU_DELTA ) )
  {
    param_map.insert(
      boost::unordered_map<std::string,std::string>::value_type(
	nnt::s_ILU_DELTA,
	zones[0].getProperty( nnt::s_ILU_DELTA )
      )
    );
  }
  else
  {
    param_map.insert(
      boost::unordered_map<std::string,std::string>::value_type(
	nnt::s_ILU_DELTA,
	"50"
      )
    );
  }

#ifdef DEBUG
  param_map.insert(
    boost::unordered_map<std::string,std::string>::value_type(
      nnt::s_ITER_SOLVER_DEBUG,
      "yes"
    )
  );
#endif

  return param_map;

}

//##############################################################################
// set_zones().
//##############################################################################

void
set_zones( std::vector<nnt::Zone>& zones, nnt::Zone& param_zone )
{

  for( size_t i = 0; i < zones.size(); i++ )
  {
    Libnucnet__Zone__iterateOptionalProperties(
      param_zone.getNucnetZone(),
      NULL,
      NULL,
      NULL,
      (Libnucnet__Zone__optional_property_iterate_function)
         nnt::copy_properties,
      zones[i].getNucnetZone()
    );
  }

}

//##############################################################################
// my_links_function().
//##############################################################################

void
my_links_function(
  const char * s_name,
  const char * s_tag1,
  const char * s_tag2,
  const char * s_value,
  link_map_t& link_map
)
{

  if( strcmp( s_name, "link" ) )
  {
    std::cerr << "Property is not a link." << std::endl;
    exit( EXIT_FAILURE );
  }

  size_t i_link = boost::lexical_cast<size_t>( s_tag1 );

  if( link_map.find( i_link )  == link_map.end() )
  {
    link_map[i_link] = boost::make_tuple( "0", "0", "0", 0. );
  }

  if( strcmp( s_tag2, "label1" ) == 0 )
    boost::get<0>( link_map.find( i_link )->second ) = s_value;
  else if( strcmp( s_tag2, "label2" ) == 0 )
    boost::get<1>( link_map.find( i_link )->second ) = s_value;
  else if( strcmp( s_tag2, "label3" ) == 0 )
    boost::get<2>( link_map.find( i_link )->second ) = s_value;
  else if( strcmp( s_tag2, "weight" ) == 0 )
    boost::get<3>( link_map.find( i_link )->second ) =
      boost::lexical_cast<double>( s_value );
  else
  {
    std::cerr << "Link contains non-valid data." << std::endl;
  }

}

//##############################################################################
// set_links().
//##############################################################################

void
set_links( user::zone_link_graph_t& g )
{

  user::zone_link_graph_t::vertex_iterator vi, vi_end;
  user::zone_link_graph_t::edge_descriptor e;
  user::vertex_multi_index vm;
  bool e_add;

  user::fill_multi_zone_vertex_hash( g, vm );

  boost::tie( vi, vi_end ) = boost::vertices( g );

  for( ; vi != vi_end; vi++ )
  {

    link_map_t my_link_map;

    Libnucnet__Zone__iterateOptionalProperties(
      g[*vi].getNucnetZone(),
      "link",
      NULL,
      NULL,
      (Libnucnet__Zone__optional_property_iterate_function) my_links_function,
      &my_link_map
    );

    for(
      link_map_t::iterator it = my_link_map.begin();
      it != my_link_map.end();
      it++
    )
    {

      boost::tie( e, e_add ) =
         boost::add_edge(
           *vi,
           user::get_vertex_from_multi_zone_hash(
             vm,
             it->second.get<0>(),
             it->second.get<1>(),
             it->second.get<2>()
           ),
           g
         );

      g[e].setWeight( it->second.get<3>() );

    }

  }

}

//##############################################################################
// get_matrix_and_vector().
//##############################################################################

std::pair<WnMatrix *, gsl_vector *>
get_matrix_and_vector(
  std::vector<nnt::Zone>& zones,
  gsl_vector * p_current
)
{

  std::pair<WnMatrix *, gsl_vector *> my_pair;
  user::zone_link_graph_t::vertex_iterator ui, ui_end;
  user::zone_link_graph_t::out_edge_iterator ei, ei_end;

  user::zone_link_graph_t g = user::make_zone_link_graph( zones );

  set_links( g );

  size_t i_species =
    Libnucnet__Nuc__getNumberOfSpecies(
      Libnucnet__Net__getNuc(
        Libnucnet__Zone__getNet( zones[0].getNucnetZone() )
      )
    );

  size_t i_offset = i_species + 1;

  size_t i_total_size = zones.size() * i_offset;

  my_pair.first = WnMatrix__new( i_total_size, i_total_size );

  my_pair.second = gsl_vector_calloc( i_total_size );

  nnt::species_list_t species_list =
    nnt::make_species_list(
      Libnucnet__Net__getNuc(
        Libnucnet__Zone__getNet( zones[0].getNucnetZone() )
      )
    );

  boost::tie( ui, ui_end ) = boost::vertices( g );

  for( ; ui != ui_end; ui++ )
  {

    size_t i_from = *ui;

    double d_mass_from =
      boost::lexical_cast<double>(
        zones[i_from].getProperty( nnt::s_ZONE_MASS )
      );

    boost::tie( ei, ei_end ) = boost::out_edges( *ui, g );

    for( ; ei != ei_end; ei++ )
    {

      user::zone_link_graph_t::vertex_descriptor v = boost::target( *ei, g );

      size_t i_to = v;

      WnMatrix__assignElement(
	my_pair.first,
	i_offset * i_from + 1,
	i_offset * i_from + 1,
	g[*ei].getWeight()
      );

      WnMatrix__assignElement(
	my_pair.first,
	i_offset * i_to + 1,
	i_offset * i_from + 1,
	-g[*ei].getWeight()
      );

      gsl_vector_set(
	my_pair.second,
	i_offset * i_from,
	gsl_vector_get(
	  my_pair.second,
	  i_offset * i_from
	) -
	d_mass_from * g[*ei].getWeight()
      );

      gsl_vector_set(
	my_pair.second,
	i_offset * i_to,
	gsl_vector_get(
	  my_pair.second,
	  i_offset * i_to
	) +
	d_mass_from * g[*ei].getWeight()
      );

      BOOST_FOREACH( nnt::Species species, species_list )
      {

	size_t i_species_offset =
	  Libnucnet__Species__getIndex( species.getNucnetSpecies() ) + 1;

	WnMatrix__assignElement(
	  my_pair.first,
	  ( i_offset * i_from ) + i_species_offset + 1,
	  ( i_offset * i_from ) + i_species_offset + 1,
	  g[*ei].getWeight()
	);

	WnMatrix__assignElement(
	  my_pair.first,
	  ( i_offset * i_to ) + i_species_offset + 1,
	  ( i_offset * i_from ) + i_species_offset + 1,
	  -g[*ei].getWeight()
	);

	gsl_vector_set(
	  my_pair.second,
	  ( i_from * i_offset ) + i_species_offset,
	  gsl_vector_get(
	    my_pair.second,
	    ( i_from * i_offset ) + i_species_offset
	  )
	  -
	  gsl_vector_get(
	    p_current,
	    ( i_from * i_offset ) + i_species_offset
	  ) * g[*ei].getWeight()
	);
	  
	gsl_vector_set(
	  my_pair.second,
	  ( i_to * i_offset ) + i_species_offset,
	  gsl_vector_get(
	    my_pair.second,
	    ( i_to * i_offset ) + i_species_offset
	  )
          +
          gsl_vector_get(
            p_current,
            ( i_from * i_offset ) + i_species_offset
          ) * g[*ei].getWeight()
	);

      }
        
    }

  }

  return my_pair;

}

//##############################################################################
// evolve_zones().
//##############################################################################

int
evolve_zones( std::vector<nnt::Zone>& zones, double dt )
{

  boost::unordered_map<std::string,std::string> param_map;

  param_map = set_iterative_solver_parameters( zones );

  user::safe_evolve_multi_zone(
    zones,
    boost::bind( get_matrix_and_vector, boost::ref( zones ), _1 ),
    boost::bind(
      user::check_multi_zone_mass_fractions, boost::ref( zones ), 1.e-8
    ),
    param_map,
    dt
  );

  return 1;

}

//##############################################################################
// update_time_and_timestep_in_zones().
//##############################################################################

void
update_time_and_timestep_in_zones(
  std::vector<nnt::Zone>& zones,
  double d_t,
  double d_dt
)
{

  for( size_t i = 0; i < zones.size(); i++ )
  {
    zones[i].updateProperty(
      nnt::s_TIME,
      boost::lexical_cast<std::string>( d_t )
    );

    zones[i].updateProperty(
      nnt::s_DTIME,
      boost::lexical_cast<std::string>( d_dt )
    );
  }

}

//##############################################################################
// main().
//##############################################################################

int
main( int argc, char **argv )
{

  std::vector<std::pair<std::string,std::string> > my_views;
  std::vector<nnt::Zone> all_zones;
  std::string s_output_file;
  double dt, dt_max, time;

  Libnucnet *p_my_nucnet, * p_param_nucnet;
  nnt::Zone param_zone;

  //============================================================================
  // Check input.
  //============================================================================

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {
    std::cout << std::endl;
    std::cout << argv[0] << " ../../data_pub/my_net.xml " <<
      "../../data_pub/multi_mass.xml " <<
      "../../data_pub/multi_mass_param.xml " <<
      "my_output.h5 " <<
      "\"[z <= 94]\" " <<
      "\"[count(reactant) = 1]\"" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc < 5 || argc > 8 )
  {
    fprintf(
      stderr,
      "\nUsage: %s xml_file zone_file param_xml hdf5_file nuc_xpath reac_xpath zone_xpath\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  xml_file = input network xml file\n\n"
    );
    fprintf(
      stderr,
      "  zone_file = input zone xml file\n\n"
    );
    fprintf(
      stderr,
      "  param_xml = input parameter xml file\n\n"
    );
    fprintf(
      stderr,
      "  hdf5_file = output hdf5 file\n\n"
    );
    fprintf(
      stderr,
      "  nuc_xpath = XPath to select nuclei (optional--required if reac_xpath present)\n\n"
    );
    fprintf(
      stderr,
      "  reac_xpath = XPath to select reactions (optional--required if zone_xpath present)\n\n"
    );
    fprintf(
      stderr,
      "  zone_xpath = XPath to select zones (optional)\n\n"
    );
    std::cout << "For an example usage, type " << std::endl << std::endl;
    std::cout << argv[0] << " --example" << std::endl << std::endl;
    return EXIT_FAILURE;
  }

  //============================================================================
  // Get data.
  //============================================================================

  p_my_nucnet = Libnucnet__new();

  if( argc == 5 )
    Libnucnet__Net__updateFromXml(
      Libnucnet__getNet( p_my_nucnet ),
      argv[1],
      "",
      ""
    );
  else if( argc == 6 )
    Libnucnet__Net__updateFromXml(
      Libnucnet__getNet( p_my_nucnet ),
      argv[1],
      argv[5],
      ""
    );
  else if( argc >= 7 )
    Libnucnet__Net__updateFromXml(
      Libnucnet__getNet( p_my_nucnet ),
      argv[1],
      argv[5],
      argv[6]
    );

  if( argc == 8 )
    Libnucnet__assignZoneDataFromXml( p_my_nucnet, argv[2], argv[7] );
  else
    Libnucnet__assignZoneDataFromXml( p_my_nucnet, argv[2], "" );

  //============================================================================
  // Sort species.
  //============================================================================

  Libnucnet__Nuc__setSpeciesCompareFunction(
    Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
    (Libnucnet__Species__compare_function) nnt::species_sort_function
  );

  Libnucnet__Nuc__sortSpecies(
    Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) )
  );

  //============================================================================
  // Get all the zones.
  //============================================================================

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  all_zones = user::get_vector_of_zones( p_my_nucnet );

  //============================================================================
  // Get parameter zone.
  //============================================================================

  p_param_nucnet = Libnucnet__new();

  Libnucnet__assignZoneDataFromXml( p_param_nucnet, argv[3], "" );

  nnt::zone_list_t param_zone_list = nnt::make_zone_list( p_param_nucnet );

  param_zone = *(param_zone_list.begin());

  //============================================================================
  // Set zones.
  //============================================================================
  
  set_zones( all_zones, param_zone );

  //============================================================================
  // Set zone views
  //============================================================================
  
  my_views.push_back( std::make_pair( "", "" ) );

  user::add_views_to_zones( all_zones, my_views );

  //============================================================================
  // Initial network limiting
  //============================================================================
  
  if( param_zone.hasProperty( nnt::s_SMALL_ABUNDANCES_THRESHOLD ) )
  {
    double d_min_abund =
      boost::lexical_cast<double>(
        param_zone.getProperty( nnt::s_SMALL_ABUNDANCES_THRESHOLD )
      );
    user::limit_zone_networks( all_zones, d_min_abund );
    user::multi_zone_zero_small_abundances( all_zones, d_min_abund );
  }
  else
  {
    user::limit_zone_networks( all_zones );
  }

  //============================================================================
  // Create the output.
  //============================================================================
  
  user::hdf5::create_output( argv[4], p_my_nucnet );

  //============================================================================
  // Evolve.
  //============================================================================
  
  if( param_zone.hasProperty( nnt::s_DTMAX ) )
  {
    dt_max =
      boost::lexical_cast<double>( param_zone.getProperty( nnt::s_DTMAX ) );
  }

  if( param_zone.hasProperty( nnt::s_TIME ) )
  {
    time =
      boost::lexical_cast<double>( param_zone.getProperty( nnt::s_TIME ) );
  }
  else
  {
    time = 0.;
  }

  dt =
    boost::lexical_cast<double>( param_zone.getProperty( nnt::s_DTIME ) );

  int i_steps = 0;

  while(
    time <
      boost::lexical_cast<double>( param_zone.getProperty( nnt::s_TEND ) )
  )
  {

    time += dt;

    update_time_and_timestep_in_zones( all_zones, time, dt );

    std::cout << "dt = " << dt << " time = " << time << " tend = " <<
      boost::lexical_cast<double>( param_zone.getProperty( nnt::s_TEND ) )
      << std::endl;

    if( !evolve_zones( all_zones, dt ) )
    {
      std::cerr << "Unable to find matrix solution." << std::endl;
      return EXIT_FAILURE;
    }

    dt = user::get_new_timestep_from_zones( all_zones, dt, 0.15, 0.15, 1.e-10 );

    if( param_zone.hasProperty( nnt::s_DTMAX ) )
    {
      if( dt > dt_max ) dt = dt_max;
    }

    if( param_zone.hasProperty( nnt::s_SMALL_ABUNDANCES_THRESHOLD ) )
    {
      double d_min_abund =
        boost::lexical_cast<double>(
          param_zone.getProperty( nnt::s_SMALL_ABUNDANCES_THRESHOLD )
      );
      user::limit_zone_networks( all_zones, d_min_abund );
      user::multi_zone_zero_small_abundances( all_zones, d_min_abund );
    }

    if(
      i_steps++ %
        boost::lexical_cast<int>( param_zone.getProperty( nnt::s_STEPS ) )
       == 0
    )
    {
      user::printout_abundances_in_zones( all_zones, time, dt );
      user::hdf5::append_zones( argv[4], p_my_nucnet );
    }
  }

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__free( p_param_nucnet );
  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
