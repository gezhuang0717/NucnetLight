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
//! \brief Example code for running a network calculation with entropy
//!        generation.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <Libnucnet.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "nnt/two_d_weak_rates.h"
#include "nnt/write_output_xml.h"

#include "user/remove_duplicate.h"
#include "user/user_rate_functions.h"
#include "user/flow_utilities.h"
#include "user/network_limiter.h"
#include "user/weak_utilities.h"
#include "user/evolve.h"
#include "user/network_utilities.h"

//##############################################################################
// Define some parameters.
//##############################################################################

#define D_REG_T        0.15    /* Time step change regulator for dt update */
#define D_REG_Y        0.15    /* Abundance change regulator for dt update */
#define D_Y_MIN_DT     1.e-10  /* Smallest y for dt update */
#define S_SOLVER       nnt::s_ARROW   
                               /* Solver type: ARROW or GSL */

//##############################################################################
// Validation.  "no" = no validation, "yes" = validation.
//##############################################################################

#define VALIDATE       "no"

//##############################################################################
// Network parameters.
//##############################################################################

#define EVOLVE_NETWORK     1
#define GENERATE_ENTROPY   1

//##############################################################################
// Strings.
//##############################################################################

#define S_RADIUS_0     "r_0"
#define S_MASS_0       "m_0"
#define S_P_0          "P_0"
#define S_RHO_0        "rho_0"

#define S_X            "x"
#define S_Y            "y"

#define S_S_NU_DOT     "neutrino entropy loss rate"
#define S_S_PHOTON_DOT "photon entropy loss rate"
#define S_T9_STEPS     "t9 steps"
#define S_T9_CUTOFF    "1.e-2"
#define S_T9_CUTOFF_FACTOR \
                       "0.1"
#define S_STEPS_BELOW  "steps below"
#define S_STEPS_ABOVE  "steps above"
#define S_ALPHA        "alpha"
#define S_FULL_RUN     "full run"
#define S_DETAILED_WEAK_RATES  "detailed weak rates"
#define S_USE_APPROXIMATE_WEAK_RATES  "use approximate weak rates"

//##############################################################################
// Prototypes.
//##############################################################################

int
func( double, const double *, double *, void * );

double
check_qse( nnt::Zone );

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  int i_step, i_steps, k = 0, i_status;
  double d_tend, d_t, d_h, d_dt_nuc;
  Libnucnet *p_my_nucnet = NULL, *p_my_output;
  Libnucnet__Reac * p_weak_lab_rates;
  Libnucnet__ReacView * p_weak_view;
  Libnucnet__Reaction * p_reaction;
  nnt::Zone zone;
  char s_property[32];
  double d_y[3];

  //============================================================================
  // Check input.
  //============================================================================

  if ( argc < 5 || argc > 7 ) {
    fprintf(
      stderr,
      "\nUsage: %s net_file zone_file out_file rho_0 xpath_nuc xpath_reac\n\n",
      argv[0]
    );
    fprintf(
      stderr, "  net_file = input network data xml filename\n\n"
    );
    fprintf(
      stderr, "  zone_file = input zone data xml filename\n\n"
    );
    fprintf(
      stderr, "  out_file = output data xml filename\n\n"
    );
    fprintf(
      stderr, "  rho_0 = initial rho\n\n"
    );
    fprintf(
      stderr,
      "  xpath_nuc = nuclear xpath expression (optional--required if xpath_reac specified)\n\n"
    );
    fprintf(
      stderr, "  xpath_reac = reaction xpath expression (optional)\n\n"
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

  if( argc == 5 ) {
    Libnucnet__Net__updateFromXml(
      Libnucnet__getNet( p_my_nucnet ),
      argv[1], 
      NULL,
      NULL
    );
  } else if( argc == 6 ) {
    Libnucnet__Net__updateFromXml(
      Libnucnet__getNet( p_my_nucnet ),
      argv[1], 
      argv[5],
      NULL
    );
  } else {
    Libnucnet__Net__updateFromXml(
      Libnucnet__getNet( p_my_nucnet ),
      argv[1], 
      argv[5],
      argv[6]
    );
  }

  Libnucnet__assignZoneDataFromXml( p_my_nucnet, argv[2], NULL );

  //============================================================================
  // Register rate functions.
  //============================================================================

  user::register_my_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

  //============================================================================
  // Set the zone.
  //============================================================================

  zone.setNucnetZone(
    Libnucnet__getZoneByLabels( p_my_nucnet, "0", "0", "0" )
  );

  //============================================================================
  // Get lab weak rates for updating at low T.
  //============================================================================

  p_weak_lab_rates =
    Libnucnet__Reac__new_from_xml(
      argv[1],
      "[single_rate]"
    );

  //============================================================================
  // Use approximate weak rates or not.
  //============================================================================

  if( zone.hasProperty( S_USE_APPROXIMATE_WEAK_RATES ) )
  {

    if( zone.getProperty( S_USE_APPROXIMATE_WEAK_RATES ) == "yes" )
    {

      user::aa522a25__update_net( Libnucnet__getNet( p_my_nucnet ) );

      p_weak_view =
	Libnucnet__ReacView__new(
	  Libnucnet__Net__getReac(
	    Libnucnet__Zone__getNet( zone.getNucnetZone() )
	  ),
	  "[source='aa522a25']"
	);

      nnt::reaction_list_t reaction_list =
	nnt::make_reaction_list(
	  Libnucnet__ReacView__getReac( p_weak_view )
	);

      BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
      {

	p_reaction =
	  Libnucnet__Reac__getReactionByString(
	    p_weak_lab_rates,
	    Libnucnet__Reaction__getString( reaction.getNucnetReaction() )
	  );

	if( p_reaction )
	{

	  Libnucnet__Reaction__updateUserRateFunctionProperty(
	    reaction.getNucnetReaction(),
	    nnt::s_LAB_RATE,
	    NULL,
	    NULL,
	    boost::lexical_cast<std::string>(
	      Libnucnet__Reaction__computeRate(
		p_reaction,
		1.,
		NULL
	      )
	    ).c_str()
	  );

        }
        else
        {

	  Libnucnet__Reaction__updateUserRateFunctionProperty(
	    reaction.getNucnetReaction(),
	    nnt::s_LAB_RATE,
	    NULL,
	    NULL,
            "0."
	  );

        }


	Libnucnet__Reaction__updateUserRateFunctionProperty(
	  reaction.getNucnetReaction(),
	  nnt::s_LAB_RATE_T9_CUTOFF,
	  NULL,
	  NULL,
	  S_T9_CUTOFF
	);
     
	Libnucnet__Reaction__updateUserRateFunctionProperty(
	  reaction.getNucnetReaction(),
	  nnt::s_LAB_RATE_T9_CUTOFF_FACTOR,
	  NULL,
	  NULL,
	  S_T9_CUTOFF_FACTOR
	);
     
      }

      Libnucnet__ReacView__free( p_weak_view );

    }
    
  }

  //============================================================================
  // Use detailed weak rates or not.  Store the lab rates for use at low T.
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

    p_weak_view =
      Libnucnet__ReacView__new(
        Libnucnet__Net__getReac(
          Libnucnet__Zone__getNet( zone.getNucnetZone() )
        ),
        "[user_rate/@key = 'two-d weak rates' or \
          user_rate/@key = 'two-d weak rates log10 ft']"
      );

    nnt::reaction_list_t reaction_list =
      nnt::make_reaction_list(
        Libnucnet__ReacView__getReac( p_weak_view )
      );

    BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
    {

      p_reaction =
        Libnucnet__Reac__getReactionByString(
          p_weak_lab_rates,
          Libnucnet__Reaction__getString( reaction.getNucnetReaction() )
        );

      if( p_reaction )
      {

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          reaction.getNucnetReaction(),
          nnt::s_LAB_RATE,
          NULL,
          NULL,
          boost::lexical_cast<std::string>(
            Libnucnet__Reaction__computeRate(
              p_reaction,
              1.,
              NULL
            )
          ).c_str()
        );

      }
      else
      {

        Libnucnet__Reaction__updateUserRateFunctionProperty(
          reaction.getNucnetReaction(),
          nnt::s_LAB_RATE,
          NULL,
          NULL,
          "0."
        );

     }

     Libnucnet__Reaction__updateUserRateFunctionProperty(
       reaction.getNucnetReaction(),
       nnt::s_LAB_RATE_T9_CUTOFF,
       NULL,
       NULL,
       S_T9_CUTOFF
     );
   
     Libnucnet__Reaction__updateUserRateFunctionProperty(
       reaction.getNucnetReaction(),
       nnt::s_LAB_RATE_T9_CUTOFF_FACTOR,
       NULL,
       NULL,
       S_T9_CUTOFF_FACTOR
     );

    }

    Libnucnet__ReacView__free( p_weak_view );

  }

  //============================================================================
  // Done with the lab rates collection, so free it.
  //============================================================================

  Libnucnet__Reac__free( p_weak_lab_rates );

  //============================================================================
  // Set two-d weak hash.
  //============================================================================

  user::set_two_d_weak_rates_hashes(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );
    
  //============================================================================
  // Set neutrino entropy loss hash if desired.
  //============================================================================

  if( zone.hasProperty( "include neutrino entropy loss" ) )
  {
    if( zone.getProperty( "include neutrino entropy loss" ) == "yes" )
      user::set_two_d_weak_energy_loss_hash(
        Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
      );
  }

  //============================================================================
  // Set the guess function.
  //============================================================================

  zone.setGuessFunction(
    (nnt::guessFunction) user::my_log10_t9_guess,
    &zone
  );

  //============================================================================
  // Remove duplicate reactions.
  //============================================================================

  user::remove_duplicate_reactions( Libnucnet__getNet( p_my_nucnet ) );

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
  // Set the initial density.
  //============================================================================

  zone.updateProperty(
    S_RHO_0,
    argv[4]
  );

  //============================================================================
  // Create output.
  //============================================================================

  p_my_output = nnt::create_output( p_my_nucnet );

  Libnucnet__setZoneCompareFunction(
    p_my_output,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  //============================================================================
  // Initialize the system.
  //============================================================================

  zone.updateProperty(
    S_RADIUS_0,
    boost::lexical_cast<std::string>(
      pow(
	(
	  3. *
	  boost::lexical_cast<double>( zone.getProperty( S_MASS_0 ) ) *
	  GSL_CONST_CGSM_SOLAR_MASS /
	  4. /
	  M_PI /
	  boost::lexical_cast<double>( zone.getProperty( S_RHO_0 ) )
	),
	1. / 3.
      )
    )
  );
    
  zone.updateProperty(
    nnt::s_RADIUS,
    zone.getProperty( S_RADIUS_0 )
  );

  d_tend =
    boost::lexical_cast<double>( zone.getProperty( nnt::s_TEND ) );

  zone.normalizeAbundances();

  d_h =
    boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) );

  if( zone.hasProperty( nnt::s_TIME ) )
    d_t =
      boost::lexical_cast<double>( zone.getProperty( nnt::s_TIME ) );
  else
    d_t = 0.0;

  zone.updateProperty(
    nnt::s_T9,
    zone.getProperty( nnt::s_T9_0 )
  );

  zone.updateProperty(
    nnt::s_RHO_0,
    argv[4]
  );

  zone.updateProperty(
    nnt::s_RHO,
    zone.getProperty( nnt::s_RHO_0 )
  );

  d_y[0] =
    boost::lexical_cast<double>( zone.getProperty( S_X ) );
  d_y[1] =
    boost::lexical_cast<double>( zone.getProperty( S_Y ) );
  d_y[2] =
    user::compute_thermo_quantity(
      zone,
      nnt::s_ENTROPY_PER_NUCLEON,
      nnt::s_TOTAL
    );

  zone.updateProperty(
    S_P_0,
    boost::lexical_cast<std::string>(
      user::compute_thermo_quantity(
        zone,
        nnt::s_PRESSURE,
        nnt::s_TOTAL
      )
    )
  );

  user::update_my_rate_functions_data( zone );  

  zone.updateProperty(
    nnt::s_DTIME,
    boost::lexical_cast<std::string>( 0. )
  );

  zone.updateProperty(
    nnt::s_YE,
    boost::lexical_cast<std::string>(
      Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 )
    )
  );

  zone.updateProperty(
    nnt::s_SDOT,
    boost::lexical_cast<std::string>(
      user::compute_entropy_change_rate(
        zone,
        Libnucnet__Zone__getEvolutionNetView( zone.getNucnetZone() )
      )
    )
  );

  const gsl_odeiv2_step_type * p_T = gsl_odeiv2_step_rkf45;

  gsl_odeiv2_step * p_step = gsl_odeiv2_step_alloc( p_T, 3 );

  gsl_odeiv2_control * p_c = gsl_odeiv2_control_y_new( 1e-4, 0.0 );

  gsl_odeiv2_evolve * p_e = gsl_odeiv2_evolve_alloc( 3 );
     
  gsl_odeiv2_system sys = {func, NULL, 3, &zone};

  user::limit_evolution_network( zone );

  //============================================================================
  // Evolve network while t < final t.
  //============================================================================

  i_step = 0;

  while(
    d_t < d_tend &&
    boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) < 20 )
  {

  //============================================================================
  // Set time.
  //============================================================================

    zone.updateProperty(
      nnt::s_TIME,
      boost::lexical_cast<std::string>( d_t )
    );  

  //============================================================================
  // Evolve step.  d_h returns with suggested next time step.
  //============================================================================

    zone.updateProperty(
      nnt::s_DT_SAFE,
      boost::lexical_cast<std::string>( d_h )
    );
  
    i_status =
      gsl_odeiv2_evolve_apply(
        p_e,
        p_c,
        p_step,
        &sys, 
        &d_t,
        d_tend,
        &d_h,
        d_y
      );
       
    if( i_status != GSL_SUCCESS )
    {
      fprintf( stderr, "error, return value=%d\n", i_status );
      return EXIT_FAILURE;
    }
       
    d_dt_nuc =
      d_t -
      boost::lexical_cast<double>( zone.getProperty( nnt::s_TIME ) );
  
  //============================================================================
  // Update properties.
  //============================================================================

    zone.updateProperty(
      nnt::s_TIME,
      boost::lexical_cast<std::string>( d_t )
    );
  
    zone.updateProperty(
      nnt::s_ENTROPY_PER_NUCLEON,
      boost::lexical_cast<std::string>( d_y[2] )
    );
  
    zone.updateProperty(
      nnt::s_RHO,
      boost::lexical_cast<std::string>(
        boost::lexical_cast<double>( zone.getProperty( S_RHO_0 ) ) *
          gsl_pow_3( 1. / d_y[0] )
      )
    );
  
  //============================================================================
  // Evolve abundances.
  //============================================================================

    if( EVOLVE_NETWORK )
    {
  
      zone.updateProperty(
        nnt::s_DTIME,
        boost::lexical_cast<std::string>( d_dt_nuc )
      );
  
      user::update_my_rate_functions_data( zone );  
  
      zone.updateProperty(
        nnt::s_T9,
        boost::lexical_cast<std::string>(
          pow(
            10.,
            zone.computeRootFromQuantity(
              (nnt::quantityFunction)
              user::my_log10_t9_from_entropy_root,
              &zone
            )
          )
  	)  
      );
  
      user::safe_evolve( zone, d_dt_nuc, d_dt_nuc );
  
      zone.updateProperty(
        nnt::s_SDOT,
        boost::lexical_cast<std::string>(
          user::compute_entropy_change_rate(
            zone,
            Libnucnet__Zone__getEvolutionNetView( zone.getNucnetZone() )
          )
        )
      );
  
      if( zone.getProperty( "include neutrino entropy loss" ) == "yes" )
      {
  
        zone.updateProperty(
          S_S_NU_DOT,
          boost::lexical_cast<std::string>(
            user::compute_approximate_neutrino_entropy_loss_rate(
              zone
            )
          )
        );
  
      }
  
    }
  
        std::cout <<
          d_t << " " <<
          zone.getProperty( nnt::s_T9 ) << " " <<
          zone.getProperty( nnt::s_RHO ) << " " <<
          d_y[0] << " " << d_y[1] << " " << d_y[2] <<
        std::endl;

    zone.updateProperty(
      S_X,
      boost::lexical_cast<std::string>( d_y[0] )
    );

    zone.updateProperty(
      S_Y,
      boost::lexical_cast<std::string>( d_y[1] )
    );

    zone.updateProperty(
      nnt::s_ENTROPY_PER_NUCLEON,
      boost::lexical_cast<std::string>( d_y[2] )
    );

    zone.updateProperty(
      nnt::s_YE,
      boost::lexical_cast<std::string>(
        Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 )
      )
    );

    zone.updateProperty(
      nnt::s_RADIUS,
      boost::lexical_cast<std::string>(
        d_y[0] *
        boost::lexical_cast<double>( zone.getProperty( S_RADIUS_0 ) )
      )
    );

  //============================================================================
  // Print out abundances.
  //============================================================================

    if(
       boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) >
       boost::lexical_cast<double>( zone.getProperty( S_T9_STEPS ) )
    )
       i_steps =
         boost::lexical_cast<int>( zone.getProperty( S_STEPS_ABOVE ) );
    else
       i_steps =
         boost::lexical_cast<int>( zone.getProperty( S_STEPS_BELOW ) );

    if( ( i_step++ % i_steps ) == 0 || d_t >= d_tend )
    {
      sprintf( s_property, "%d", ++k );
      Libnucnet__relabelZone(
        p_my_nucnet,
        zone.getNucnetZone(),
        s_property,
        NULL,
        NULL
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
      &d_dt_nuc,
      D_REG_T,
      D_REG_Y,
      D_Y_MIN_DT
    );

    if( d_dt_nuc < d_h ) d_h = d_dt_nuc;

    if ( d_t + d_h > d_tend ) {

      d_h = d_tend - d_t;

    }

    zone.normalizeAbundances();

    user::limit_evolution_network( zone );

  }  

  //============================================================================
  // Write output.
  //============================================================================

  Libnucnet__updateZoneXmlMassFractionFormat(
    p_my_output,
    "%.15e"
  );

  Libnucnet__writeToXmlFile( p_my_output, argv[3] );

  //============================================================================
  // Clean up and exit.
  //============================================================================

  gsl_odeiv2_evolve_free( p_e );
  gsl_odeiv2_control_free( p_c );
  gsl_odeiv2_step_free( p_step );

  Libnucnet__free( p_my_nucnet );
  Libnucnet__free( p_my_output );

  return EXIT_SUCCESS;

}

//##############################################################################
// func().
//##############################################################################

int
func( double d_t, const double d_y[], double d_f[], void * p_params )
{

  double d_dt, d_dt_fac;
  gsl_vector * p_abundances, * p_abundance_changes;

  nnt::Zone zone = *( nnt::Zone * ) p_params;

  p_abundances = Libnucnet__Zone__getAbundances( zone.getNucnetZone() );

  p_abundance_changes =
    Libnucnet__Zone__getAbundanceChanges( zone.getNucnetZone() );

  d_dt =
    d_t 
    -
    boost::lexical_cast<double>( zone.getProperty( nnt::s_TIME ) );

  d_dt_fac =
    fabs(
      d_dt /
      boost::lexical_cast<double>( zone.getProperty( nnt::s_DT_SAFE ) )
    );

  zone.updateProperty(
    nnt::s_DTIME,
    boost::lexical_cast<std::string>( d_dt )
  );

  zone.updateProperty(
    nnt::s_ENTROPY_PER_NUCLEON,
    boost::lexical_cast<std::string>( d_y[2] )
  );

  zone.updateProperty(
    nnt::s_RHO,
    boost::lexical_cast<std::string>(
      boost::lexical_cast<double>(
        zone.getProperty(
          S_RHO_0
        )
      ) * gsl_pow_3( 1. / d_y[0] )
    )
  );

  if( EVOLVE_NETWORK && d_dt_fac > 1.e-10 )
  {

    user::update_my_rate_functions_data( zone );  

    if( zone.getProperty( S_FULL_RUN ) == "yes" )
    {
      zone.updateProperty(
	nnt::s_T9,
        boost::lexical_cast<std::string>(
	  pow(
	    10.,
	    zone.computeRootFromQuantity(
	      (nnt::quantityFunction)
		user::my_log10_t9_from_entropy_root,
	      &zone
	    )  
	  )
        )
      );
    }
    else
    {
      zone.updateProperty(
	nnt::s_T9,
        boost::lexical_cast<std::string>(
	  pow(
	    10.,
	    zone.computeRootFromQuantity(
	      (nnt::quantityFunction)
		user::compute_log10_t9_entropy_root,
	      &zone
	    )  
	  )
        )
      );
    }

    user::safe_evolve( zone, d_dt, d_dt );

  }

  if( GENERATE_ENTROPY )
  {

    if( d_dt_fac > 1.e-10 )
    {

      user::update_my_rate_functions_data( zone );  

      d_f[2] =
        user::compute_entropy_change_rate(
          zone,
          Libnucnet__Zone__getEvolutionNetView( zone.getNucnetZone() )
        );

   }
    else
    {

      d_f[2] =
        boost::lexical_cast<double>(
          zone.getProperty( nnt::s_SDOT )
        );

    }

    if( zone.getProperty( "include neutrino entropy loss" ) == "yes" )
    {

      zone.updateProperty(
        S_S_NU_DOT,
        boost::lexical_cast<std::string>(
          user::compute_approximate_neutrino_entropy_loss_rate(
            zone
          )
        )
      );

      d_f[2] -=
        boost::lexical_cast<double>(
          zone.getProperty( S_S_NU_DOT )
        );

         std::cout << "s_nu_dot = " << zone.getProperty( S_S_NU_DOT ) << "; ";
    }

    if( zone.getProperty( "include photon entropy loss" ) == "yes" )
    {

      zone.updateProperty(
        S_S_PHOTON_DOT,
        boost::lexical_cast<std::string>(
          user::compute_photon_entropy_loss_rate( zone )
        )
      );

      d_f[2] -=
        boost::lexical_cast<double>(
          zone.getProperty( S_S_PHOTON_DOT )
        );

         std::cout << "s_photon_dot = " << zone.getProperty( S_S_PHOTON_DOT );
         std::cout << std::endl;
    }

  }
  else
    d_f[2] = 0;

  d_f[0] = d_y[1];

  d_f[1] = 
    boost::lexical_cast<double>( zone.getProperty( S_ALPHA ) ) *
    (
      boost::lexical_cast<double>( zone.getProperty( S_P_0 ) ) /
      (
        boost::lexical_cast<double>( zone.getProperty( S_RHO_0 ) )
        *
        gsl_pow_2( 
          boost::lexical_cast<double>( zone.getProperty( S_RADIUS_0 ) )
        )
      )
    )
    *
    (
      (
        user::compute_thermo_quantity(
          zone,
	  nnt::s_PRESSURE,
	  nnt::s_TOTAL
        ) /
        boost::lexical_cast<double>( zone.getProperty( S_P_0 ) )
      ) * gsl_pow_2( d_y[0] )
      -
      1. / gsl_pow_2( d_y[0] )
    );

   printf(
     "%s  %s  %s  %g  %g\n",
     zone.getProperty( nnt::s_DTIME ).c_str(),
     zone.getProperty( nnt::s_T9 ).c_str(),
     zone.getProperty( nnt::s_RHO ).c_str(),
     d_y[2],
     d_f[2]
   );

  Libnucnet__Zone__updateAbundances( zone.getNucnetZone(), p_abundances );
      
  Libnucnet__Zone__updateAbundanceChanges(
    zone.getNucnetZone(),
    p_abundance_changes
  );

  gsl_vector_free( p_abundances );
  gsl_vector_free( p_abundance_changes );
      
  return GSL_SUCCESS;

} 

//##############################################################################
// check_qse().
//##############################################################################

double
check_qse(
  nnt::Zone zone
)
{

  Libnuceq * p_equil;
  Libnuceq__Cluster * p_cluster;
  double d_result;

  if(
    boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) )
    <
    4
  )
    return 1.e6;

  p_equil =
    Libnuceq__new(
      Libnucnet__Net__getNuc( Libnucnet__Zone__getNet( zone.getNucnetZone() ) )
    );

  Libnuceq__setYe(
    p_equil,
    Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 )
  );

  p_cluster =
    Libnuceq__newCluster( p_equil, "[z >= 6]" );

  Libnuceq__Cluster__updateConstraint(
    p_cluster,
    user::compute_cluster_abundance_moment(
      zone,
      "[z >= 6]",
      "a",
      0
    )
  );
   
  Libnuceq__computeEquilibrium(
    p_equil,
    boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ),
    boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) )
  );

  d_result = Libnuceq__Cluster__getMukT( p_cluster );

  Libnuceq__free( p_equil );

  return d_result;

}
