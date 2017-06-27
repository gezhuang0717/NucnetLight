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
//! \brief Example code for running a multiple-zone network calculation
//!          with OpenMP.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <omp.h>
#include <vector>
#include <Libnucnet.h>

#include "user/evolve.h"
#include "user/network_utilities.h"

#include "nnt/string_defs.h"
#include "nnt/auxiliary.h"
#include "nnt/iter.h"
#include "nnt/write_output_xml.h"

#include "user/remove_duplicate.h"
#include "user/screen.h"
#include "user/coul_corr.h"
#include "user/user_rate_functions.h"

//##############################################################################
// Define some parameters.
//##############################################################################

#define D_DT0          1.e-15  /* Initial timestep */

//##############################################################################
// Validation.  "no" = no validation, "yes" = validation.
//##############################################################################

#define VALIDATE       "no"

//##############################################################################
// Strings.
//##############################################################################

#define S_SOLVER_TYPE   nnt::s_ARROW  /* nnt::s_ARROW or nnt::s_GSL */
#define S_USE_APPROXIMATE_WEAK_RATES  "use approximate weak rates"
#define S_PARAMS "params"
#define S_NEW_RATES_XML "new rates xml"

//##############################################################################
// Screening.  "no" = no screening, "yes" = screening.
//##############################################################################

#define S_USE_SCREENING  "yes"

//##############################################################################
// Other defines.
//##############################################################################

#define I_TAUS  100

//##############################################################################
// set_zones().
//##############################################################################

void
set_zones(
  Libnucnet * p_nucnet,
  std::vector<double>& tau_vector,
  std::vector<nnt::Zone>& zone_vector
)
{

  nnt::Zone param_zone, zone;
  user::view_multi views;

  //============================================================================
  // Get the parameter zone.
  //============================================================================

  param_zone.setNucnetZone(
    Libnucnet__getZoneByLabels( p_nucnet, S_PARAMS, "0", "0" )
  );

  //============================================================================
  // Update reactions, if desired.
  //============================================================================

  if( param_zone.hasProperty( S_NEW_RATES_XML ) )
  {
    Libnucnet__Reac__updateFromXml(
      Libnucnet__Net__getReac( Libnucnet__getNet( p_nucnet ) ),
      param_zone.getProperty( S_NEW_RATES_XML ).c_str(),
      NULL
    );
  }

  //============================================================================
  // Add rate modification views.
  //============================================================================

  Libnucnet__Zone__iterateOptionalProperties(
    param_zone.getNucnetZone(),
    NULL,
    nnt::s_RATE_MODIFICATION_VIEW,
    NULL,
    (Libnucnet__Zone__optional_property_iterate_function)
       user::modify_rates_views,
    &views
  );

  BOOST_FOREACH( user::view view, views )
  {
    Libnucnet__Zone__updateNetView(
      param_zone.getNucnetZone(),
      view.nuc_xpath.c_str(),
      view.reac_xpath.c_str(),
      NULL,
      Libnucnet__NetView__copy(
        param_zone.getNetView(
          view.nuc_xpath.c_str(),
          view.reac_xpath.c_str()
        )
      )
    );
  }
 
  //============================================================================
  // Set the other net views.
  //============================================================================

  param_zone.getNetView(
    "",
    ""
  );

  if( param_zone.hasProperty( nnt::s_LOW_T_NETWORK_REAC_XPATH ) )
  {
    Libnucnet__Zone__updateNetView(
      param_zone.getNucnetZone(),
      "",
      param_zone.getProperty( nnt::s_LOW_T_NETWORK_REAC_XPATH ).c_str(),
      NULL,
      Libnucnet__NetView__new( 
        Libnucnet__Zone__getNet( param_zone.getNucnetZone() ),
        "",
        param_zone.getProperty( nnt::s_LOW_T_NETWORK_REAC_XPATH ).c_str()
      )
    );
  }

  //============================================================================
  // Set the solver.
  //============================================================================

  param_zone.updateProperty(
    nnt::s_SOLVER,
    S_SOLVER_TYPE
  );

  param_zone.updateProperty(
    nnt::s_ARROW_WIDTH,
    "3"
  );

  //============================================================================
  // Create and add zones.  Notice that, since the parameter zone is already
  // present in p_my_nucnet, we swap its label before adding the new zone,
  // then we relabel the new zone, and swap back the label of the parameter
  // zone.
  //============================================================================

  int i = 0;

  BOOST_FOREACH( double tau, tau_vector )
  {

    zone.setNucnetZone(
      Libnucnet__Zone__copy( param_zone.getNucnetZone() )
    );

    Libnucnet__relabelZone(
      p_nucnet,
      param_zone.getNucnetZone(),
      "save",
      NULL,
      NULL
    );

    if( !Libnucnet__addZone( p_nucnet, zone.getNucnetZone() ) )
    {
      std::cerr << "Couldn't add zone with tau = " << tau << std::endl;
      exit( EXIT_FAILURE );
    }

    Libnucnet__relabelZone(
      p_nucnet,
      zone.getNucnetZone(),
      boost::lexical_cast<std::string>( i++ ).c_str(),
      NULL,
      NULL
    );

    Libnucnet__relabelZone(
      p_nucnet,
      param_zone.getNucnetZone(),
      S_PARAMS,
      NULL,
      NULL
    );

    zone.updateProperty(
      nnt::s_TAU,
      boost::lexical_cast<std::string>( tau )
    );

    zone.updateProperty(
      nnt::s_T9,
      param_zone.getProperty( nnt::s_T9_0 )
    );

    zone.updateProperty(
      nnt::s_RHO,
      param_zone.getProperty( nnt::s_RHO_0 )
    );

    if( !zone.hasProperty( nnt::s_TIME ) )
      zone.updateProperty(
        nnt::s_TIME,
        "0."
      );

    if( !zone.hasProperty( nnt::s_DTIME ) )
      zone.updateProperty(
        nnt::s_DTIME,
        boost::lexical_cast<std::string>( D_DT0 )
      );

    zone_vector.push_back( zone );
	  
    user::limit_evolution_network( zone );

  }

}

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  size_t i;
  std::vector<double> tau_vector;
  std::vector<nnt::Zone> zone_vector;
  Libnucnet__Zone * p_zone;
  Libnucnet *p_my_nucnet;

  //============================================================================
  // Check input.
  //============================================================================

  if ( argc != 4 && argc != 5 )
  {
    fprintf(
      stderr,
      "\nUsage: %s net_xml zone_xml out_xml nuc_xpath\n\n",
      argv[0]
    );
    fprintf(
      stderr, "  net_xml = network data xml filename\n\n"
    );
    fprintf(
      stderr, "  zone_xml = zone data xml filename\n\n"
    );
    fprintf(
      stderr, "  out_xml = output xml filename\n\n"
    );
    fprintf(
      stderr, "  nuc_xpath = XPath expression for nuclei\n\n"
    );
    return EXIT_FAILURE;
  }

  //============================================================================
  // Validate input file.
  //============================================================================

  if( strcmp( VALIDATE, "yes" ) == 0 )
  {
    if( !Libnucnet__is_valid_input_xml( argv[1] ) ) {
      fprintf( stderr, "Not valid libnucnet input!\n" );
      return EXIT_FAILURE;
    }
  }

  //============================================================================
  // Read and store input.
  //============================================================================

  p_my_nucnet = Libnucnet__new();

  if( argc == 5 )
    Libnucnet__Net__updateFromXml(
      Libnucnet__getNet( p_my_nucnet ),
      argv[1],
      argv[4],
      NULL
    );  
  else
    Libnucnet__Net__updateFromXml(
      Libnucnet__getNet( p_my_nucnet ),
      argv[1],
      NULL,
      NULL
    );  

  Libnucnet__assignZoneDataFromXml(
    p_my_nucnet,
    argv[2],
    NULL
  );

  //============================================================================
  // Register user-supplied rate functions.
  //============================================================================

  user::register_my_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

  //============================================================================
  // Remove duplicate reactions.
  //============================================================================

  user::remove_duplicate_reactions( Libnucnet__getNet( p_my_nucnet ) );

  //============================================================================
  // Sort the nuclei if using the arrow solver.
  //============================================================================

  if(
    strcmp( S_SOLVER_TYPE, nnt::s_ARROW ) == 0
  )
  {

    Libnucnet__Nuc__setSpeciesCompareFunction(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
      (Libnucnet__Species__compare_function) nnt::species_sort_function
    );

    Libnucnet__Nuc__sortSpecies(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ) 
    );

  }

  //============================================================================
  // Get the taus.
  //============================================================================
  
  for( int i = 1; i <= I_TAUS; i++ )
  {

     tau_vector.push_back(
       static_cast<double>( i ) / 100.
     );

  }

  //============================================================================
  // Set the zones.
  //============================================================================
  
  p_zone = Libnucnet__getZoneByLabels( p_my_nucnet, "0", "0", "0" );

  Libnucnet__relabelZone(
    p_my_nucnet,
    p_zone,
    S_PARAMS,
    NULL,
    NULL
  );

  set_zones( p_my_nucnet, tau_vector, zone_vector );

  Libnucnet__removeZone( p_my_nucnet, p_zone );

  //============================================================================
  // Evolve zones.
  //============================================================================

  #pragma omp parallel for schedule( dynamic, 1 )
    for(
      i = 0;
      i < zone_vector.size();
      i++
    )
    {

      user::evolve_zone(
        zone_vector[i],
        boost::lexical_cast<double>( zone_vector[i].getProperty( nnt::s_TEND ) )
      );

      Libnucnet__Zone__clearRates( zone_vector[i].getNucnetZone() );

      Libnucnet__Zone__clearNetViews( zone_vector[i].getNucnetZone() );

    }
  
  //============================================================================
  // Write out.
  //============================================================================

  Libnucnet__writeToXmlFile( p_my_nucnet, argv[3] );

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
