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
//! \brief Example code for running a single zone network calculation.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <Libnucnet.h>

#include "nnt/two_d_weak_rates.h"
#include "nnt/write_output_xml.h"

#include "user/remove_duplicate.h"
#include "user/user_rate_functions.h"
#include "user/network_limiter.h"
#include "user/network_utilities.h"
#include "user/rate_modifiers.h"
#include "user/evolve.h"
#include "user/hydro.h"
#include "user/hdf5_routines.h"

//##############################################################################
// Define some parameters.
//##############################################################################

#define D_DT0          1.e-20       // Initial time step
#define D_REG_T        0.15         // Time step change regulator for dt update
#define D_REG_Y        0.15         // Abundance change regulator for dt update 
#define D_Y_MIN_DT     1.e-10       // Smallest y for dt update
#define S_SOLVER       nnt::s_ARROW // Solver type: ARROW or GSL

#define S_DETAILED_WEAK_RATES  "detailed weak rates"
#define S_FLOW_CURRENT_XML_FILE  "flow current xml file"
#define S_INTEGRATED_CURRENTS  "integrated currents"

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  int i_step;
  double d_t, d_dt;
  Libnucnet *p_my_nucnet = NULL, *p_flow_current_nucnet = NULL;
  nnt::Zone zone, flow_current_zone;
  std::set<std::string> isolated_species_set;

  //============================================================================
  // Get the nucnet.
  //============================================================================

  p_my_nucnet = user::get_nucnet( argc, argv );

  //============================================================================
  // Register rate functions.
  //============================================================================

  user::register_my_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

  //============================================================================
  // Get the zone.
  //============================================================================

  if( !user::set_zone( p_my_nucnet, zone, argv ) )
  {
    std::cerr << "Couldn't set zone."  << std::endl;
    return EXIT_FAILURE;
  }

  //============================================================================
  // Initialize time.
  //============================================================================

  if( zone.hasProperty( nnt::s_DTIME ) )
    d_dt = boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) );
  else
  {
    d_dt = D_DT0;
    zone.updateProperty(
      nnt::s_DTIME,
      boost::lexical_cast<std::string>( d_dt )
    );
  }

  if( zone.hasProperty( nnt::s_TIME ) )
    d_t = boost::lexical_cast<double>( zone.getProperty( nnt::s_TIME ) );
  else
  {
    d_t = 0;

    zone.updateProperty(
      nnt::s_TIME,
      boost::lexical_cast<std::string>( d_t )
    );

  }

  //============================================================================
  // Initialize zone.
  //============================================================================

  user::initialize_zone( zone, argv );

  //============================================================================
  // Use approximate weak rates or not.
  //============================================================================

  if(
    zone.hasProperty( nnt::s_USE_APPROXIMATE_WEAK_RATES ) &&
    zone.getProperty( nnt::s_USE_APPROXIMATE_WEAK_RATES ) == "yes"
  )
    user::aa522a25__update_net( Libnucnet__getNet( p_my_nucnet ) );

  //============================================================================
  // Update with detailed weak rates.
  //============================================================================

  if( zone.hasProperty( S_DETAILED_WEAK_RATES ) )
  {

    Libnucnet__Reac__updateFromXml(
      Libnucnet__Net__getReac(
        Libnucnet__Zone__getNet( zone.getNucnetZone() )
      ),
      zone.getProperty( S_DETAILED_WEAK_RATES ).c_str(),
      NULL
    );

    user::set_two_d_weak_rates_hashes(
      Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
    );

  }

  //============================================================================
  // If not set, set the neutrino chemical potential to -inf.
  //============================================================================
  
  if( !zone.hasProperty( nnt::s_MU_NUE_KT ) )
    zone.updateProperty(
      nnt::s_MU_NUE_KT,
      "-inf"
    ); 

  //============================================================================
  // Remove duplicate reactions.
  //============================================================================

  user::remove_duplicate_reactions( Libnucnet__getNet( p_my_nucnet ) );

  //============================================================================
  // Remove isolated species if desired.
  //============================================================================

  if( zone.hasProperty( nnt::s_SPECIES_REMOVAL_REAC_XPATH ) )
  {

    isolated_species_set =
      user::get_isolated_species(
        Libnucnet__getNet( p_my_nucnet ),
        "",
        zone.getProperty( nnt::s_SPECIES_REMOVAL_REAC_XPATH )
      );

    BOOST_FOREACH( std::string s_species, isolated_species_set )
    {

//    Careful that you don't remove a species with non-zero abundance!

      std::cout << s_species << std::endl;

      Libnucnet__Nuc__removeSpecies(
        Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc(
            Libnucnet__getNet( p_my_nucnet )
          ),
          s_species.c_str()
        )
      );

    }
    
  }

  //============================================================================
  // Sort the nuclei if using the arrow solver.
  //============================================================================

  if( S_SOLVER == nnt::s_ARROW )
  {

    Libnucnet__Nuc__setSpeciesCompareFunction(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
      (Libnucnet__Species__compare_function) nnt::species_sort_function
    );

    Libnucnet__Nuc__sortSpecies(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ) 
    );

    zone.updateProperty( nnt::s_SOLVER, nnt::s_ARROW );

    zone.updateProperty( nnt::s_ARROW_WIDTH, "3" );

  }

  //============================================================================
  // Create the output hdf5 file.  Do this after nuclide sort.
  //============================================================================

  user::hdf5::create_output( 
    zone.getProperty( S_OUTPUT_FILE ).c_str(),
    p_my_nucnet
  );

  //============================================================================
  // Create current nucnet.
  //============================================================================

  if( zone.hasProperty( S_FLOW_CURRENT_XML_FILE ) )
  {

    p_flow_current_nucnet = nnt::create_output( p_my_nucnet );

    flow_current_zone.setNucnetZone(
      Libnucnet__Zone__new(
        Libnucnet__getNet( p_flow_current_nucnet ),
        S_INTEGRATED_CURRENTS,
        NULL,
        NULL
      )
    );

    Libnucnet__addZone(
      p_flow_current_nucnet,
      flow_current_zone.getNucnetZone()
    );

    user::copy_zone_abundances_as_properties(
      zone,
      flow_current_zone,
      nnt::s_INITIAL_ABUNDANCE
    );

  }

  //============================================================================
  // Print out the rate modification views.
  //============================================================================

  user::print_modified_reactions( zone );

  //============================================================================
  // Limit evolution network.
  //============================================================================

  user::limit_evolution_network( zone );

  //============================================================================
  // Initialize exposure for specific species.
  //============================================================================

  if( zone.hasProperty( nnt::s_SPECIFIC_SPECIES ) )
  {
    zone.updateProperty(
      nnt::s_EXPOSURE,
      zone.getProperty( nnt::s_SPECIFIC_SPECIES ),
      "0"
    );
  }

  //============================================================================
  // Evolve network while t < final t.
  //============================================================================

  i_step = 0;

  while ( d_t < boost::lexical_cast<double>( zone.getProperty( nnt::s_TEND ) ) )
  {

    d_t += d_dt;

  //============================================================================
  // Set dt and t.
  //============================================================================

    zone.updateProperty(
      nnt::s_DTIME,
      boost::lexical_cast<std::string>( d_dt )
    );

    zone.updateProperty(
      nnt::s_TIME,
      boost::lexical_cast<std::string>( d_t )
    );

  //============================================================================
  // Update temperature and density.  Update dt and time again in case
  // changed in update_zone_properties.
  //============================================================================

    user::update_zone_properties( zone );

    d_t = boost::lexical_cast<double>( zone.getProperty( nnt::s_TIME ) );

    d_dt = boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) );

  //============================================================================
  // Evolve abundances.
  //============================================================================

    user::evolve( zone );

  //============================================================================
  // Update exposure and, if desired, network currents.
  //============================================================================

    user::update_exposures( zone );

    if( zone.hasProperty( S_FLOW_CURRENT_XML_FILE ) )
    {
      user::update_flow_currents( zone, flow_current_zone );
    }

  //============================================================================
  // Print out abundances.
  //============================================================================

    if(
       (
         i_step % 
           boost::lexical_cast<int>( zone.getProperty( nnt::s_STEPS ) ) == 0 ||
         d_t >= boost::lexical_cast<double>( zone.getProperty( nnt::s_TEND ) )
       )
    )
    {
      user::hdf5::append_zones(
        zone.getProperty( S_OUTPUT_FILE ).c_str(),
        p_my_nucnet
      );
      zone.printAbundances();
    }

  //============================================================================
  // Update timestep.
  //============================================================================

    Libnucnet__Zone__updateTimeStep(
      zone.getNucnetZone(),
      &d_dt,
      D_REG_T,
      D_REG_Y,
      D_Y_MIN_DT
    );

    if( boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) > 10. )
      zone.normalizeAbundances();

    if(
      d_t + d_dt >
        boost::lexical_cast<double>( zone.getProperty( nnt::s_TEND ) )
    )
    {
      d_dt =
        boost::lexical_cast<double>( zone.getProperty( nnt::s_TEND ) ) - d_t;
    }

    user::limit_evolution_network( zone );

    i_step++;

  }  

  //============================================================================
  // Write output.
  //============================================================================

  if( zone.hasProperty( S_FLOW_CURRENT_XML_FILE ) )
  {

    user::copy_zone_abundances_as_properties(
      zone,
      flow_current_zone,
      nnt::s_FINAL_ABUNDANCE
    );

    Libnucnet__writeToXmlFile(
      p_flow_current_nucnet,
      zone.getProperty( S_FLOW_CURRENT_XML_FILE ).c_str()
    );

    Libnucnet__free( p_flow_current_nucnet );

  }
        
  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
