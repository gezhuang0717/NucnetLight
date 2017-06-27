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
//! \brief Multi-zone example code.
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
// get_matrix_and_vector().
//##############################################################################

std::pair<WnMatrix *, gsl_vector *>
get_matrix_and_vector(
  std::vector<nnt::Zone>& zones,
  gsl_vector * p_current
)
{

  std::pair<WnMatrix *, gsl_vector *> my_pair;

  //============================================================================
  // Get species list.
  //============================================================================
  
  nnt::species_list_t species_list =
    nnt::make_species_list(
      Libnucnet__Net__getNuc(
        Libnucnet__Zone__getNet( zones[0].getNucnetZone() )
      )
    );

  size_t i_species = species_list.size();

  size_t i_offset = i_species;

  size_t i_total_size = zones.size() * i_offset;

  my_pair.first = WnMatrix__new( i_total_size, i_total_size );

  my_pair.second = gsl_vector_calloc( i_total_size );

  for( size_t i_from = 0; i_from < zones.size(); i_from++ )
  {

    std::vector<size_t> i_to_vector;
   
    if( i_from > 0 ) i_to_vector.push_back( i_from - 1 );
    if( i_from < zones.size() - 1 ) i_to_vector.push_back( i_from + 1 );
   
    BOOST_FOREACH( size_t i_to, i_to_vector )
    {

      BOOST_FOREACH( nnt::Species species, species_list )
      {

	size_t i_species_offset =
	  Libnucnet__Species__getIndex( species.getNucnetSpecies() );

	WnMatrix__assignElement(
	  my_pair.first,
	  ( i_offset * i_from ) + i_species_offset + 1,
	  ( i_offset * i_from ) + i_species_offset + 1,
	  D_MIX_RATE
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
	  ) * D_MIX_RATE
	);
	  
	WnMatrix__assignElement(
	  my_pair.first,
	  ( i_offset * i_from ) + i_species_offset + 1,
	  ( i_offset * i_to ) + i_species_offset + 1,
	  -D_MIX_RATE
	);

	gsl_vector_set(
	  my_pair.second,
	  ( i_from * i_offset ) + i_species_offset,
	  gsl_vector_get(
	    my_pair.second,
	    ( i_from * i_offset ) + i_species_offset
	  )
	  +
	  gsl_vector_get(
	    p_current,
	    ( i_to * i_offset ) + i_species_offset
	  ) * D_MIX_RATE
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

  Libnucnet *p_my_nucnet, * p_param_nucnet;
  nnt::Zone param_zone;

  //============================================================================
  // Check input.
  //============================================================================

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {
    std::cout << std::endl;
    std::cout << argv[0] << " ../../data_pub/my_net.xml " <<
      "../../data_pub/multi_zone.xml " <<
      "../../data_pub/multi_zone_param.xml " <<
      "my_output.h5 " <<
      "\"[z <= 30]\"" <<
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
  
  double dt =
    boost::lexical_cast<double>( param_zone.getProperty( nnt::s_DTIME ) );
  double d_time = 0.;

  int i_steps = 0;

  while(
    d_time <
      boost::lexical_cast<double>( param_zone.getProperty( nnt::s_TEND ) )
  )
  {

    d_time += dt;

    update_time_and_timestep_in_zones( all_zones, d_time, dt );

    std::cout << "dt = " << dt << " time = " << d_time << " tend = " <<
      boost::lexical_cast<double>( param_zone.getProperty( nnt::s_TEND ) )
      << std::endl;

    if( !evolve_zones( all_zones, dt ) )
    {
      std::cerr << "Unable to find matrix solution." << std::endl;
      return EXIT_FAILURE;
    }

    dt = user::get_new_timestep_from_zones( all_zones, dt, 0.15, 0.15, 1.e-10 );

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
      user::printout_abundances_in_zones( all_zones, d_time, dt );
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
