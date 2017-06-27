//////////////////////////////////////////////////////////////////////////////
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
//! \brief Example code for running a network calculation at constant entropy.
////////////////////////////////////////////////////////////////////////////////

/*##############################################################################
// Includes.
//############################################################################*/

#include <Libnucnet.h>
#include "nnt/auxiliary.h"
#include "nnt/string_defs.h"
#include "nnt/write_output_xml.h"
#include "user/evolve.h"

#include "user/user_rate_functions.h"
#include "user/neutrino_rate_functions.h"
#include "user/network_utilities.h"
#include "user/remove_duplicate.h"

#include "my_user/my_hydro.h"
/*##############################################################################
// Define some parameters.
//############################################################################*/

#define MY_BUF_SIZE    32      /* String buffer size */
#define D_DT0          1.e-05  /* Initial time step */
#define D_REG_T        0.15    /* Time step change regulator for dt update */
#define D_REG_Y        0.15    /* Abundance change regulator for dt update */
#define D_Y_MIN_DT     1.e-10  /* Smallest y for dt update */

/*##############################################################################
// Validation.  "no" = no validation, "yes" = validation.
//############################################################################*/

#define VALIDATE       "no"

/*##############################################################################
// Strings.
//############################################################################*/

#define S_SOLVER_TYPE   nnt::s_ARROW  /* nnt::s_ARROW or nnt::s_GSL */
#define S_USE_APPROXIMATE_WEAK_RATES  "use approximate weak rates"

/*##############################################################################
// main().
//############################################################################*/

int main( int argc, char * argv[] ) {

  int i_step;
  double d_tend, d_t, d_dt;
  Libnucnet *p_my_nucnet, *p_my_output;

  nnt::Zone zone;

  //============================================================================
  // Check input.
  //============================================================================

  if ( argc < 4 || argc > 7 ) {
    fprintf(
      stderr,
      "\nUsage: %s net_file zone_xml out_xml xpath_nuc xpath_reac neutrino_xml\n\n",
      argv[0]
    );
    fprintf(
      stderr, "  net_file = nuclear network data xml filename\n\n"
    );
    fprintf(
      stderr, "  zone_file = input zone data xml filename\n\n"
    );
    fprintf(
      stderr, "  out_xml = output xml filename\n\n"
    );
    fprintf(
      stderr,
      "  xpath_nuc = nuclear xpath expression (optional--required if xpath_reac specified)\n\n"
    );
    fprintf(
      stderr, "  xpath_reac = reaction xpath expression (optional--required if neutrino_xml specified)\n\n"
    );
    fprintf(
      stderr, "  neutrino_xml = neutrino input xml filename (optional)\n\n"
    );
    return EXIT_FAILURE;
  }

  //============================================================================
  // Validate input file.
  //============================================================================

  if( strcmp( VALIDATE, "yes" ) == 0 )
  {
    if( !Libnucnet__Net__is_valid_input_xml( argv[1] ) ) {
      fprintf( stderr, "Not valid libnucnet input!\n" );
      return EXIT_FAILURE;
    }
  }

  //============================================================================
  // Read and store input.
  //============================================================================

  p_my_nucnet = Libnucnet__new();

  if( argc == 4 )
    Libnucnet__Net__updateFromXml(
      Libnucnet__getNet( p_my_nucnet ),
      argv[1],
      NULL,
      NULL
    );  
  else if( argc == 5 )
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
      argv[4],
      argv[5]
    );  

  //============================================================================
  // Input zone data.
  //============================================================================

  Libnucnet__assignZoneDataFromXml(
    p_my_nucnet,
    argv[2],
    NULL
  );

  //============================================================================
  // Get the zone.
  //============================================================================

  zone.setNucnetZone(
    Libnucnet__getZoneByLabels( p_my_nucnet, "0", "0", "0" )
  );

  if( !zone.getNucnetZone() )
  {
    fprintf( stderr, "No input zone with labels (0,0,0).\n" );
    exit( EXIT_FAILURE );
  }


  //============================================================================
  // Update with neutrino data.
  //============================================================================

  if( argc == 7 )
  {

    Libnucnet__Reac__updateFromXml(
      Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) ),
      argv[6],
      NULL
    );

    user::set_nu_nucl_hash(
      Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
    );

  }

  //============================================================================
  // Register user-supplied rate functions.
  //============================================================================

  user::register_my_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

  user::register_my_neutrino_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

  //============================================================================
  // Remove duplicate reactions.
  //============================================================================

  user::remove_duplicate_reactions( Libnucnet__getNet( p_my_nucnet ) );

  //============================================================================
  // Update with approximate weak rates if desired.
  //============================================================================

  if( zone.hasProperty( S_USE_APPROXIMATE_WEAK_RATES ) )
  {
    if( zone.getProperty( S_USE_APPROXIMATE_WEAK_RATES ) == "yes" )
      user::aa522a25__update_net( Libnucnet__getNet( p_my_nucnet ) );
  }

  //============================================================================
  // If not set, set the neutrino chemical potential to -inf.
  //============================================================================
  
  if( !zone.hasProperty( nnt::s_MU_NUE_KT ) )
  {
    zone.updateProperty(
      nnt::s_MU_NUE_KT,
      "-inf"
    );
  }

  //============================================================================
  // Create output network structure.
  //============================================================================

  p_my_output = nnt::create_output( p_my_nucnet );

  //============================================================================
  // Set solver.
  //============================================================================

  Libnucnet__Zone__updateProperty(
    zone.getNucnetZone(), nnt::s_SOLVER, NULL, NULL, S_SOLVER_TYPE
  );

  //============================================================================
  // Set the guess function.
  //============================================================================

  zone.setGuessFunction(
    (nnt::guessFunction) user::my_log10_t9_guess,
    &zone
  );

  //============================================================================
  // Sort the nuclei if using the arrow solver.
  //============================================================================

  if( zone.getProperty( nnt::s_SOLVER ) == nnt::s_ARROW )
  {

    Libnucnet__Nuc__setSpeciesCompareFunction(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
      (Libnucnet__Species__compare_function)
         nnt::species_sort_function
    );

    Libnucnet__Nuc__sortSpecies(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ) 
    );

    Libnucnet__Nuc__setSpeciesCompareFunction(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_output ) ),
      (Libnucnet__Species__compare_function)
         nnt::species_sort_function
    );

    Libnucnet__Nuc__sortSpecies(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_output ) )
    );

    zone.updateProperty( nnt::s_ARROW_WIDTH, "3" );

  }

  //============================================================================
  // Start with network in NSE at input t9 and Ye.
  //============================================================================

  zone.updateProperty(
    nnt::s_RHO,
    boost::lexical_cast<std::string>(
      pow(
        10.,
        zone.computeRootFromQuantity(
          (nnt::quantityFunction)
            user::compute_log10_density_entropy_root_with_equilibrium,
          &zone
        )
      )
    )
  );

  zone.updateProperty(
    nnt::s_RHO_0,
    zone.getProperty( nnt::s_RHO )
  );

  //============================================================================
  // Limit network.
  //============================================================================

  user::limit_evolution_network( zone );

  //============================================================================
  // Evolve network while t < final t.
  //============================================================================

  i_step = 0;
  d_t =
    boost::lexical_cast<double>(
      zone.getProperty(
        nnt::s_TIME
      )
    );
  d_tend =
    boost::lexical_cast<double>(
      zone.getProperty(
        nnt::s_TEND
      )
    );

  while( d_t < d_tend )
  {

    d_dt =
      boost::lexical_cast<double>(
        zone.getProperty( nnt::s_DTIME )
      );

    d_t += d_dt;

    zone.updateProperty(
      nnt::s_TIME,
      boost::lexical_cast<std::string>( d_t )
    );

  //============================================================================
  // Update density and radius.  Swap neutrinos if below resonance.
  //============================================================================

    my_user::update_zone_properties( zone );

    if( zone.hasProperty( S_RHO_RES ) )
      user::swap_neutrinos( zone );

  //============================================================================
  // Evolve abundances.
  //============================================================================

    zone.updateProperty(
      nnt::s_T9,
      boost::lexical_cast<std::string>(
        pow(
          10.,
          zone.computeRootFromQuantity(
            (nnt::quantityFunction) user::my_log10_t9_from_entropy_root,
            &zone
          )
        )
      )
    );

    user::evolve( zone );

  //============================================================================
  // Print out abundances.
  //============================================================================

    if(
      (
        i_step %
        boost::lexical_cast<int>(
          zone.getProperty(
            nnt::s_STEPS
          )
        ) == 0 
      ) ||
        d_t >= d_tend
    )
    {
      zone.updateProperty(
        nnt::s_YE,
        boost::lexical_cast<std::string>(
          Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 )
        )
      );
      zone.printAbundances();
      nnt::write_xml( p_my_output, zone.getNucnetZone() );
      Libnucnet__writeToXmlFile( p_my_output, argv[3] );
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

    if ( d_t + d_dt > d_tend ) {

      d_dt = d_tend - d_t;

    }

    zone.updateProperty(
      nnt::s_DTIME,
      boost::lexical_cast<std::string>( d_dt )
    );

    user::limit_evolution_network( zone );

    zone.normalizeAbundances();

    i_step++;

  }  

  //============================================================================
  // Write out and free p_my_output.
  //============================================================================

  Libnucnet__setZoneCompareFunction(
    p_my_output,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );
  Libnucnet__writeToXmlFile( p_my_output, argv[3] );
  Libnucnet__free( p_my_output );

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__free( p_my_nucnet );
  return EXIT_SUCCESS;

}
