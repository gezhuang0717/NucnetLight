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
//! \brief Example code for running a network calculation with energy
//!        generation.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <Libnucnet.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include <boost/format.hpp>

#include "nnt/two_d_weak_rates.h"
#include "nnt/write_output_xml.h"
#include "user/remove_duplicate.h"
#include "user/user_rate_functions.h"
#include "user/network_limiter.h"
#include "user/flow_utilities.h"

#include "user/evolve.h"

//##############################################################################
// Define some parameters.
//##############################################################################

#define D_DT0          1.e-15  /* Initial time step */
#define D_REG_T        0.15    /* Time step change regulator for dt update */
#define D_REG_Y        0.15    /* Abundance change regulator for dt update */
#define D_Y_MIN_DT     1.e-10  /* Smallest y for dt update */
#define I_SOLVER       nnt::s_ARROW   
                               /* Solver type: ARROW or GSL */
#define S_USE_APPROXIMATE_WEAK_RATES "use approximate weak rates"

#define D_ALPHA        1.

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
#define S_P_0          "P_0"
#define S_RHO_0        "rho_0"

#define S_X            "x"
#define S_Y            "y"

#define S_XDOT         "xdot"
#define S_YDOT         "ydot"

#define S_PARTICLE     nnt::s_TOTAL
#define S_ENERGY_GENERATION   "energy generation"
#define S_COMPUTE_ENERGY  "compute energy"

#define S_REAC_XPATH   "[                      \
                          reactant = 'c12' or  \
                          reactant = 'o16' or  \
                          reactant = 'ne20' or \
                          reactant = 'mg24' or \
                          reactant = 'si28' or \
                          reactant = 'p31'  or \
                          (reactant = 'he4' and product = 'c12') \
                       ]" 

//##############################################################################
// Prototypes.
//##############################################################################

int
func( double, const double *, double *, void * );

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  int i_step, k = 0, i_status;
  double d_tend, d_t, d_h, d_dt_nuc;
  Libnucnet *p_my_nucnet = NULL, *p_my_output;
  nnt::Zone zone;
  char s_property[32];
  double d_y[3];

  //============================================================================
  // Check input.
  //============================================================================

  if ( argc < 9 || argc > 11 ) {
    fprintf(
      stderr,
      "\nUsage: %s in_file out_file t9 rho radius t_end m_step mu_nue_kT xpath_nuc xpath_reac\n\n",
      argv[0]
    );
    fprintf(
      stderr, "  in_file = input single zone data xml filename\n\n"
    );
    fprintf(
      stderr, "  out_file = output data xml filename\n\n"
    );
    fprintf(
      stderr, "  t9 = input temperature\n\n"
    );
    fprintf(
      stderr, "  rho = input density (in g/cc)\n\n"
    );
    fprintf(
      stderr, "  radius = initial radius\n\n"
    );
    fprintf(
      stderr, "  t_end = duration to evolve system\n\n"
    );
    fprintf(
      stderr, "  m_step = frequency of steps to print out\n\n"
    );
    fprintf(
      stderr, "  mu_nue_kT = neutrino chemical potential / kT\n\n"
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

  if( argc == 9 ) {
    p_my_nucnet = Libnucnet__new_from_xml( argv[1], NULL, NULL, NULL );
  } else if( argc == 10 ) {
    p_my_nucnet = Libnucnet__new_from_xml( argv[1], argv[9], NULL, NULL );
  } else {
    p_my_nucnet = Libnucnet__new_from_xml( argv[1], argv[9], argv[10], NULL );
  }

  d_tend = atof( argv[6] );

  //============================================================================
  // Register rate functions.
  //============================================================================

  user::register_my_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

  //============================================================================
  // Set the zone and normalize the abundances.
  //============================================================================

  zone.setNucnetZone(
    Libnucnet__getZoneByLabels( p_my_nucnet, "0", "0", "0" )
  );

  zone.normalizeAbundances();

  //============================================================================
  // Use approximate weak rates or not.
  //============================================================================

  if( zone.hasProperty( S_USE_APPROXIMATE_WEAK_RATES ) )
  {
    if( zone.getProperty( S_USE_APPROXIMATE_WEAK_RATES ) == "yes" )
    {
      user::aa522a25__update_net( Libnucnet__getNet( p_my_nucnet ) );
    }
  }

  //============================================================================
  // Set the neutrino chemical potential.
  //============================================================================
  
  Libnucnet__Zone__updateProperty(
    zone.getNucnetZone(),
    nnt::s_MU_NUE_KT,
    NULL,
    NULL,
    argv[8]
  ); 

  //============================================================================
  // Remove duplicate reactions.
  //============================================================================

  user::remove_duplicate_reactions( Libnucnet__getNet( p_my_nucnet ) );

  //============================================================================
  // Sort the nuclei if using the arrow solver.
  //============================================================================

  if( I_SOLVER == nnt::s_ARROW )
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
  // Create output.
  //============================================================================

  p_my_output = nnt::create_output( p_my_nucnet );

  Libnucnet__setZoneCompareFunction(
    p_my_output,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  //============================================================================
  // Set initial temperature and density.
  //============================================================================

  Libnucnet__Zone__updateProperty(
    zone.getNucnetZone(),
    nnt::s_T9_0,
    NULL,
    NULL,
    argv[3]
  );
  
  Libnucnet__Zone__updateProperty(
    zone.getNucnetZone(),
    nnt::s_T9,
    NULL,
    NULL,
    argv[3]
  );
  
  Libnucnet__Zone__updateProperty(
    zone.getNucnetZone(),
    nnt::s_RHO_0,
    NULL,
    NULL,
    argv[4]
  );
  
  Libnucnet__Zone__updateProperty(
    zone.getNucnetZone(),
    nnt::s_RHO,
    NULL,
    NULL,
    argv[4]
  );
  
  //============================================================================
  // Initialize the system.
  //============================================================================

  d_h = D_DT0;
  d_t = 0.0;

  d_y[0] = 1.;
  d_y[1] = 0.;
  d_y[2] = atof( argv[3] );

  zone.updateProperty(
    S_P_0,
    boost::lexical_cast<std::string>(
      user::compute_thermo_quantity(
        zone,
        nnt::s_PRESSURE,
        S_PARTICLE
      )
    )
  );

  user::update_my_rate_functions_data( zone );

  zone.updateProperty(
    S_ENERGY_GENERATION,
    boost::lexical_cast<std::string>(
      user::compute_energy_generation_rate_per_nucleon(
        zone,
        Libnucnet__Zone__getEvolutionNetView( zone.getNucnetZone() )
      )
    )
  );

  Libnucnet__Zone__updateProperty(
    zone.getNucnetZone(),
    S_RADIUS_0,
    NULL,
    NULL,
    argv[5]
  );

  Libnucnet__Zone__updateProperty(
    zone.getNucnetZone(),
    S_COMPUTE_ENERGY,
    NULL,
    NULL,
    "no"
  );

  const gsl_odeiv2_step_type * p_T = gsl_odeiv2_step_rk8pd;

  gsl_odeiv2_step * p_step = gsl_odeiv2_step_alloc( p_T, 3 );

  gsl_odeiv2_control * p_c = gsl_odeiv2_control_y_new( 1e-4, 0.0 );

  gsl_odeiv2_evolve * p_e = gsl_odeiv2_evolve_alloc( 3 );
     
  gsl_odeiv2_system sys = {func, NULL, 3, &zone};

  user::limit_evolution_network( zone );

  //============================================================================
  // Evolve network while t < final t.
  //============================================================================

  i_step = 0;

  while ( d_t < d_tend )
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
    d_t - boost::lexical_cast<double>( zone.getProperty( nnt::s_TIME ) );

  //============================================================================
  // Update properties.
  //============================================================================

    zone.updateProperty(
      nnt::s_TIME,
      boost::lexical_cast<std::string>( d_t )
    );

    zone.updateProperty(
      nnt::s_T9,
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

    user::safe_evolve( zone, d_dt_nuc, d_dt_nuc );

  }

        std::cout <<
          d_t << " " <<
          zone.getProperty( nnt::s_T9 ) << " " <<
          d_y[0] << " " << d_y[1] << " " << d_y[2] << " " <<
          user::compute_thermo_quantity(
            zone,
            nnt::s_ENTROPY_PER_NUCLEON,
            S_PARTICLE
          ) <<
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
      S_ENERGY_GENERATION,
      boost::lexical_cast<std::string>(
        user::compute_energy_generation_rate_per_nucleon(
          zone,
          Libnucnet__Zone__getEvolutionNetView( zone.getNucnetZone() )
        )
      )
    );

    Libnucnet__Zone__updateProperty(
      zone.getNucnetZone(),
      S_COMPUTE_ENERGY,
      NULL,
      NULL,
      "no"
    );

  //============================================================================
  // Print out abundances.
  //============================================================================

    if( ( i_step++ % atoi( argv[7] ) ) == 0 || d_t >= d_tend )
    {
      sprintf( s_property, "%d", ++k );
      Libnucnet__relabelZone(
        p_my_nucnet,
        zone.getNucnetZone(),
        s_property,
        NULL,
        NULL
      );
      zone.printAbundances( );
      nnt::write_xml( p_my_output, zone.getNucnetZone() );
      Libnucnet__writeToXmlFile( p_my_output, argv[2] );
      user::limit_evolution_network( zone );
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

//    zone.normalizeAbundances();

  }  

  //============================================================================
  // Write output.
  //============================================================================

  Libnucnet__updateZoneXmlMassFractionFormat(
    p_my_output,
    "%.15e"
  );

  Libnucnet__writeToXmlFile( p_my_output, argv[2] );

  //============================================================================
  // Clean up and exit.
  //============================================================================

  gsl_odeiv2_evolve_free( p_e );
  gsl_odeiv2_control_free( p_c );
  gsl_odeiv2_step_free( p_step );

  Libnucnet__free( p_my_nucnet );
  return EXIT_SUCCESS;

}

//##############################################################################
// func().
//##############################################################################

int
func( double d_t, const double d_y[], double d_f[], void * p_params )
{

  double d_dt, d_energy_generation, d_ds;
  gsl_vector * p_abundances, * p_abundance_changes;

  nnt::Zone zone = *( nnt::Zone * ) p_params;

  p_abundances = Libnucnet__Zone__getAbundances( zone.getNucnetZone() );

  p_abundance_changes =
    Libnucnet__Zone__getAbundanceChanges( zone.getNucnetZone() );

  d_dt =
    d_t 
    -
    boost::lexical_cast<double>( zone.getProperty( nnt::s_TIME ) )
    + 1.e-25;

  zone.updateProperty(
    nnt::s_DTIME,
    boost::lexical_cast<std::string>( d_dt )
  );

  zone.updateProperty(
    nnt::s_T9,
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

  if( GENERATE_ENTROPY && EVOLVE_NETWORK )
  {

    user::update_my_rate_functions_data( zone );

    if(
      zone.getProperty( S_COMPUTE_ENERGY ) == "yes"
    )
    {
      d_energy_generation =
        user::compute_energy_generation_rate_per_nucleon(
          zone,
          Libnucnet__Zone__getEvolutionNetView( zone.getNucnetZone() )
	);
    }
    else
    {
    
      user::safe_evolve( zone, d_dt, d_dt );

      d_energy_generation =
	boost::lexical_cast<double>( zone.getProperty( S_ENERGY_GENERATION ) );

      Libnucnet__Zone__updateProperty(
	zone.getNucnetZone(),
	S_COMPUTE_ENERGY,
	NULL,
	NULL,
	"yes"
      );
    }

    d_ds =
      user::compute_entropy_change_rate(
        zone,
        Libnucnet__Zone__getEvolutionNetView( zone.getNucnetZone() )
      );

  }
  else
  {
    d_energy_generation = 0.;
    d_ds = 0;
  }

  d_f[0] = d_y[1];

  d_f[1] = 
    D_ALPHA *
    (
      boost::lexical_cast<double>(
        zone.getProperty( S_P_0 )
      ) /
      (
        boost::lexical_cast<double>( zone.getProperty( S_RHO_0 ) )
        *
        gsl_pow_2( 
          boost::lexical_cast<double>(
            zone.getProperty( S_RADIUS_0 )
          )
        )
      )
    )
    *
    (
      (
        user::compute_thermo_quantity(
          zone,
	  nnt::s_PRESSURE,
          S_PARTICLE
        ) /
        boost::lexical_cast<double>(
          zone.getProperty(
            S_P_0
          )
        )
      ) * gsl_pow_2( d_y[0] )
      -
      1. / gsl_pow_2( d_y[0] )
    );

   d_f[2] =
     ( 
       1. / 
         user::compute_thermo_quantity(
           zone,
           nnt::s_SPECIFIC_HEAT_PER_NUCLEON, 
           S_PARTICLE 
         ) 
     )
     *
     (
       -boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) *
          GSL_CONST_NUM_GIGA
         * user::compute_thermo_quantity( zone, nnt::s_DPDT, S_PARTICLE )
         * 3. * gsl_pow_2( d_y[0] ) * d_y[1]
       / 
       (
         boost::lexical_cast<double>(
           zone.getProperty(
             S_RHO_0
           )
         ) *
         GSL_CONST_NUM_AVOGADRO
       )
       + d_energy_generation
     );

   d_f[2] /= GSL_CONST_NUM_GIGA;

   std::cout <<
     boost::format( "%8s %8s %8s %.5e %.5e %5e\n" ) %
     zone.getProperty( nnt::s_DTIME ) %
     zone.getProperty( nnt::s_T9 ) %
     zone.getProperty( nnt::s_RHO ) %
     d_energy_generation %
     d_ds;

   zone.updateProperty(
     S_XDOT,
     boost::lexical_cast<std::string>( d_f[0] )
   );

   zone.updateProperty(
     S_YDOT,
     boost::lexical_cast<std::string>( d_f[1] )
   );

  Libnucnet__Zone__updateAbundances( zone.getNucnetZone(), p_abundances );
      
  Libnucnet__Zone__updateAbundanceChanges(
    zone.getNucnetZone(),
    p_abundance_changes
  );
      
  return GSL_SUCCESS;

} 

