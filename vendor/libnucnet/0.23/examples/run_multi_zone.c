/*//////////////////////////////////////////////////////////////////////////////
// <file type="public">
//   <license>
//     This file was originally written by Bradley S. Meyer.
//
//     This is free software; you can redistribute it and/or modify it
//     under the terms of the GNU General Public License as published by
//     the Free Software Foundation; either version 2 of the License, or
//     (at your option) any later version.
//
//     This software is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     Please see the src/README.txt file in this distribution for more
//     information.
//   </license>
//   <description>
//     <abstract>
//       Example to demonstrate how to use libnucnet routines to create
//       a new Libnucnet structure from an input xml file, evolve a multi-zone
//       nuclear network system, and clear the structure and free the
//       allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

/*##############################################################################
// Includes.
//############################################################################*/

#include <Libnucnet.h>
#include "remove_duplicate.h"
#include "screen.h"
#include "coul_corr.h"
#include "user_rate_functions.h"

/*##############################################################################
// Define some parameters.
//############################################################################*/

#define D_DT0          1.e-10  /* Initial time step */
#define D_MIN          1.e-4   /* Newton-Raphson convergence criterion */
#define D_Y_MIN        1.e-30  /* Smallest y for convergence */
#define D_Y_MIN_PRINT  1.e-30  /* Smallest y to print out */

#define D_REGT         0.15    /* Max fractional change in dt */
#define D_REGY         0.15    /* Max fractional change in y */
#define D_YMIN         1.e-10  /* Smallest y for change in dt */

#define I_ITMAX        4       /* Maximum number of Newton-Raphson iterations */

#define T9             "t9"
#define RHO            "rho"

/*##############################################################################
// Validation.  "no" = no validation, "yes" = validation.
//############################################################################*/

#define VALIDATE       "yes"

/*##############################################################################
// Screening.  "no" = no screening, "yes" = screening.
//############################################################################*/

#define USE_SCREENING  "yes"

/*##############################################################################
// User structures.
//############################################################################*/

typedef struct {
  Libnucnet *pLibnucnet;
  gsl_vector *pYOld;
  gsl_vector *pRhs;
  gsl_vector *pSolution;
  WnMatrix *pMatrix;
  double dMixTime;
  double dDt;
  double dCheck;
  size_t iIter;
} user_data;

  user_screening_data my_screening_data;
  user_coul_corr_data my_corr_data = {-0.9052, 0.6322};

/*##############################################################################
// Prototypes.
//############################################################################*/

int evolve_system( Libnucnet *, double, double );
int print_zone( Libnucnet__Zone *, void * );
void print_abunds( Libnucnet__Zone * );
int print_abund( Libnucnet__Species *, Libnucnet__Zone * );
int get_next_timestep( Libnucnet__Zone *, double * );
int set_zone_abundances( Libnucnet__Zone *, gsl_vector * );
int set_zone( Libnucnet__Zone *, user_data * );
int set_mixing_terms( Libnucnet__Zone *, user_data * );
int update_abundances( Libnucnet__Zone *, user_data * );
int update_abundance_changes( Libnucnet__Zone *, gsl_vector * );
int initialize_zone( Libnucnet__Zone *, void * );
int check_zone_for_t9_and_rho( Libnucnet__Zone *, int * );
int set_zone_for_mixing_only( Libnucnet__Zone *, void * );
int zone_compare( const Libnucnet__Zone *, const Libnucnet__Zone * );
int add_copy_of_net_view( Libnucnet__Zone *, Libnucnet__NetView * );

/*##############################################################################
// main()
//############################################################################*/

int main( int argc, char * argv[] ) {

  double d_tend, d_t, d_dt;
  double d_mix_time;
  int i_step;
  int i_all_zones_have_t9_and_rho = 1;
  Libnucnet *p_my_nucnet;
  Libnucnet__NetView *p_view;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if ( argc < 5 || argc > 7 ) {
    fprintf(
      stderr,
      "\nUsage: %s network_file mix_time t_end mix_flag xpath_nuc xpath_reac\n\n",
      argv[0]
    );
    fprintf(
      stderr, "  network_file = input nuclear data xml filename\n\n"
    );
    fprintf(
      stderr, "  mix_time = mixing time between zones (0 = no mixing)\n\n"
    );
    fprintf(
      stderr, "  t_end = duration to evolve system\n\n"
    );
    fprintf(
      stderr, "  mix_flag = flag for full calculation (non-zero) or only mixing (0)\n\n"
    );
    fprintf(
      stderr, "  xpath_nuc = nuclear xpath expression (optional--required if xpath_reac specified\n\n"
    );
    fprintf(
      stderr, "  xpath_reac = reaction xpath expression (optional)\n\n"
    );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Validate input xml.
  //==========================================================================*/

  if( strcmp( VALIDATE, "yes" ) == 0 )
  {
    if( !Libnucnet__is_valid_input_xml( argv[1] ) ) {
      fprintf( stderr, "Not valid libnucnet input xml!\n" );
      return EXIT_FAILURE;
    }
  }

  /*============================================================================
  // Read and store network data.
  //==========================================================================*/

  if( argc == 5 ) {
    p_my_nucnet = Libnucnet__new_from_xml( argv[1], NULL, NULL, NULL );
  } else if( argc == 6 ) {
    p_my_nucnet = Libnucnet__new_from_xml( argv[1], argv[5], NULL, NULL );
  } else {
    p_my_nucnet = Libnucnet__new_from_xml( argv[1], argv[5], argv[6], NULL );
  }
 
  /*============================================================================
  // Register user-supplied rate functions.
  //==========================================================================*/

  register_my_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

  /*============================================================================
  // Check all zones for t9 and rho.
  //==========================================================================*/

  Libnucnet__iterateZones(
    p_my_nucnet,
    (Libnucnet__Zone__iterateFunction) check_zone_for_t9_and_rho,
    &i_all_zones_have_t9_and_rho
  );

  if( !i_all_zones_have_t9_and_rho )
  {
    fprintf( stderr, "Not all zones have t9 and rho!\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // If mixing only desired, set all zone temperatures to zero.
  //==========================================================================*/

  if( atoi( argv[4] ) == 0 )
    Libnucnet__iterateZones(
      p_my_nucnet,
      (Libnucnet__Zone__iterateFunction) set_zone_for_mixing_only,
      NULL
    );

  /*============================================================================
  // Iterate over zones to remove duplicates and set screening.
  //==========================================================================*/

  Libnucnet__iterateZones(
    p_my_nucnet,
    (Libnucnet__Zone__iterateFunction) initialize_zone,
    NULL
  );

  /*============================================================================
  // Iterate over zones to set valid net view.
  //==========================================================================*/

  p_view =
    Libnucnet__NetView__new( 
      Libnucnet__getNet( p_my_nucnet ),
      "",
      ""
    );

  Libnucnet__iterateZones(
    p_my_nucnet,
    (Libnucnet__Zone__iterateFunction) add_copy_of_net_view,
    p_view
  );

  Libnucnet__NetView__free( p_view );

  /*============================================================================
  // Store mixing time and final time.
  //==========================================================================*/

  d_mix_time = atof( argv[2] );
  d_tend = atof( argv[3] );

  /*============================================================================
  // Initialize time.
  //==========================================================================*/

  d_dt = D_DT0;
  d_t = 0.0;

  /*============================================================================
  // Evolve network while t < final t.
  //==========================================================================*/

  i_step = 0;

  while ( d_t < d_tend ) {

    d_t += d_dt;

    if ( evolve_system( p_my_nucnet, d_mix_time, d_dt ) != 0 ) {

      fprintf( stderr, "Error evolving system!\n" );
      return EXIT_FAILURE;

    }

  /*============================================================================
  // Print out abundances every 25th time step.  Before the zone iteration,
  // set the zone compare function to zone_compare, which compares zones
  // by the integer value of the first zone label.  After the zone iteration,
  // reset the zone compare function to the default.
  //==========================================================================*/

    if( i_step % 25 == 0 ) {

      printf( "t = %10.4e, dt = %10.4e\n\n", d_t, d_dt );

      Libnucnet__setZoneCompareFunction(
        p_my_nucnet,
        (Libnucnet__Zone__compare_function) zone_compare
      );

      Libnucnet__iterateZones(
        p_my_nucnet,
        (Libnucnet__Zone__iterateFunction) print_zone,
        NULL
      );

      Libnucnet__clearZoneCompareFunction( p_my_nucnet );

    }

  /*============================================================================
  // Update timestep.
  //==========================================================================*/

    Libnucnet__iterateZones(
      p_my_nucnet,
      (Libnucnet__Zone__iterateFunction) get_next_timestep,
      &d_dt
    ); 

    if ( d_t + d_dt > d_tend ) {

      d_dt = d_tend - d_t;

    }

    i_step++;

  }  

  /*============================================================================
  // Print out final abundances.
  //==========================================================================*/

  printf( "t = %10.4e, dt = %10.4e\n\n", d_t, d_dt );

  Libnucnet__iterateZones(
    p_my_nucnet,
    (Libnucnet__Zone__iterateFunction) print_zone,
    NULL
  );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  Libnucnet__free( p_my_nucnet );
  return EXIT_SUCCESS;

}

/*##############################################################################
// evolve_system()
//############################################################################*/

int
evolve_system( 
  Libnucnet *self, double d_mix_time, double d_dt 
) {

  size_t i_species, i_zones, i_dim, i_iter;
  gsl_vector *p_abundance_changes;
  user_data *p_user_data;
  
  /*============================================================================
  // Allocate Memory.
  //==========================================================================*/

  i_species =
    Libnucnet__Nuc__getNumberOfSpecies(
      Libnucnet__Net__getNuc(
        Libnucnet__getNet( self )
      )
    );

  i_zones = Libnucnet__getNumberOfZones( self );

  i_dim = i_species * i_zones;

  p_abundance_changes = gsl_vector_calloc( i_dim );

  p_user_data = ( user_data * ) malloc( sizeof( user_data ) );

  p_user_data->pRhs = gsl_vector_calloc( i_dim );

  p_user_data->dMixTime = d_mix_time;
  p_user_data->dDt = d_dt;
  p_user_data->pYOld = gsl_vector_calloc( i_dim );

  /*============================================================================
  // Store initial abundances.
  //==========================================================================*/

  p_user_data->pLibnucnet = self;

  Libnucnet__iterateZones(
    self,
    (Libnucnet__Zone__iterateFunction) set_zone_abundances,
    p_user_data->pYOld
  );

  /*============================================================================
  // Newton-Raphson iterations.
  //==========================================================================*/

  for( i_iter = 0; i_iter < I_ITMAX; i_iter++ ) {

    p_user_data->pMatrix = WnMatrix__new( i_dim, i_dim );

    p_user_data->iIter = i_iter;

    /*==========================================================================
    // Loop Over Zones.
    //========================================================================*/

    Libnucnet__iterateZones(
      self,
      (Libnucnet__Zone__iterateFunction) set_zone,
      p_user_data
    ); 


    /*==========================================================================
    // Assign mixing terms.
    //========================================================================*/

    if( d_mix_time < 0. ) {
       fprintf( stderr, "Invalid mixing time" );
       exit( EXIT_FAILURE );
    }

    if( d_mix_time > 0. )
      Libnucnet__iterateZones(
        self,
        (Libnucnet__Zone__iterateFunction) set_mixing_terms,
        p_user_data
      );

    /*========================================================================
    // Add 1/dt to diagonals.
    //========================================================================*/

    WnMatrix__addValueToDiagonals( p_user_data->pMatrix, 1.0 / d_dt );

    /*========================================================================
    // Solve matrix.
    //========================================================================*/

    p_user_data->pSolution =
      WnMatrix__solve( p_user_data->pMatrix, p_user_data->pRhs );
  
    /*==========================================================================
    // Update abundances.
    //========================================================================*/

    p_user_data->dCheck = 0.;

    Libnucnet__iterateZones(
      self,
      (Libnucnet__Zone__iterateFunction) update_abundances,
      p_user_data
    );

    /*==========================================================================
    // Update abundance changes.
    //========================================================================*/

    gsl_vector_add( p_abundance_changes, p_user_data->pSolution );

    /*==========================================================================
    // Free matrix and solution vector.
    //========================================================================*/

    gsl_vector_free( p_user_data->pSolution );
    WnMatrix__free( p_user_data->pMatrix );

    /*==========================================================================
    // Break if solution converged.
    //========================================================================*/

    if( p_user_data->dCheck < D_MIN ) break;

  }

  /*============================================================================
  // Update abundance changes.
  //==========================================================================*/

  Libnucnet__iterateZones(
    self,
    (Libnucnet__Zone__iterateFunction) update_abundance_changes,
    p_abundance_changes
  );


  /*==========================================================================
  // Free allocated memory.
  //========================================================================*/

  gsl_vector_free( p_abundance_changes );

  gsl_vector_free( p_user_data->pRhs );
  gsl_vector_free( p_user_data->pYOld );

  free( p_user_data );

  return 0;

}

/*##############################################################################
// print_zone()
//############################################################################*/

int
print_zone( Libnucnet__Zone *p_zone, void *p_data )
{

  if( p_data ) {
    fprintf( stderr, "Routine should have no extra data!\n" );
    return 0;
  }

  fprintf(
    stdout,
    "Zone: %s %s %s\n\n",
    Libnucnet__Zone__getLabel( p_zone, 1 ),
    Libnucnet__Zone__getLabel( p_zone, 2 ),
    Libnucnet__Zone__getLabel( p_zone, 3 )
  );

  print_abunds( p_zone );

  return 1;

}

/*##############################################################################
// print_abunds()
//############################################################################*/

void
print_abunds( Libnucnet__Zone *p_zone ) {

  Libnucnet__Nuc__iterateSpecies(
    Libnucnet__Net__getNuc(
      Libnucnet__Zone__getNet( p_zone )
    ),
    (Libnucnet__Species__iterateFunction) print_abund,
    p_zone
  );

  fprintf(
    stdout,
    "\n1 - xsum = %g\n",
    1. - Libnucnet__Zone__computeAMoment( p_zone, 1 )
  );

}

/*##############################################################################
// print_abund()
//############################################################################*/

int
print_abund( Libnucnet__Species *p_species, Libnucnet__Zone *p_zone )
{

  double d_abund;

  d_abund = Libnucnet__Zone__getSpeciesAbundance( p_zone, p_species );

  if ( d_abund > D_Y_MIN_PRINT ) {

    printf(
      "%5u%5u%14.4e%14.4e\n",
      Libnucnet__Species__getZ( p_species ),
      Libnucnet__Species__getA( p_species ),
      d_abund,
      Libnucnet__Zone__getSpeciesAbundanceChange( p_zone, p_species )
    );

  }

  return 1;

}

/*##############################################################################
// get_next_timestep().
//############################################################################*/

int
get_next_timestep( Libnucnet__Zone *p_zone, double *p_dt )
{

  Libnucnet__Zone__updateTimeStep(
    p_zone,
    p_dt,
    D_REGT,
    D_REGY,
    D_YMIN
  );

  return 1;

}

/*##############################################################################
// set_zone_abundances().
//############################################################################*/

int
set_zone_abundances( Libnucnet__Zone *p_zone, gsl_vector *p_old )
{

  gsl_vector *p_abundances;
  gsl_vector_view v_old;
  size_t i_zone;

  i_zone = (size_t) atol( Libnucnet__Zone__getLabel( p_zone, 1 ) );

  p_abundances =
    Libnucnet__Zone__getAbundances( p_zone );

  v_old =
    gsl_vector_subvector(
      p_old,
      i_zone * WnMatrix__get_gsl_vector_size( p_abundances ),
      WnMatrix__get_gsl_vector_size( p_abundances )
    );

  gsl_vector_memcpy( &v_old.vector, p_abundances );

  gsl_vector_free( p_abundances );

  return 1;

}

/*##############################################################################
// set_zone()
//############################################################################*/

int
set_zone( Libnucnet__Zone *p_zone, user_data *p_user_data )
{

  size_t i_zone;
  gsl_vector *p_abund, *p_rhs;
  gsl_vector_view v_old, v_rhs;
  WnMatrix *p_matrix;
  double d_rhoe;

  /*========================================================================
  // Get zone.
  //======================================================================*/

  i_zone = (size_t) atol( Libnucnet__Zone__getLabel( p_zone, 1 ) );

  /*========================================================================
  // If T9 <= 0, return.
  //======================================================================*/

  if(
    atof( Libnucnet__Zone__getProperty( p_zone, T9, NULL, NULL ) ) <= 0.
  )
    return 1;

  /*========================================================================
  // Update data for user-rate functions.
  //======================================================================*/

  d_rhoe =
    atof(
      Libnucnet__Zone__getProperty( p_zone, RHO, NULL, NULL )
    ) *
    Libnucnet__Zone__computeZMoment( p_zone, 1 );

  Libnucnet__Zone__updateDataForUserRateFunction(
    p_zone,
    CF88_WEAK_FIT,
    &d_rhoe
  );
      
  /*============================================================================
  // Set screening, if desired.
  //==========================================================================*/

  if( strcmp( USE_SCREENING, "yes" ) == 0 ) {

    /*==========================================================================
    // Set screening data and function.
    //========================================================================*/

    my_screening_data.dYe2 =
      Libnucnet__Zone__computeZMoment( p_zone, 2 );

    Libnucnet__Zone__setScreeningFunction(
      p_zone,
      (Libnucnet__Net__screening_function) my_screening_function,
      &my_screening_data
    );

    /*==========================================================================
    // Set correction function.
    //========================================================================*/

    Libnucnet__Zone__setNseCorrectionFactorFunction(
      p_zone,
      (Libnucnet__Species__nseCorrectionFactorFunction) my_coulomb_correction,
      &my_corr_data
    );

  }

  /*========================================================================
  // Compute rates.  
  //======================================================================*/

  Libnucnet__Zone__computeRates(
    p_zone,
    atof(
      Libnucnet__Zone__getProperty( p_zone, T9, NULL, NULL )
    ),
    atof(
      Libnucnet__Zone__getProperty( p_zone, RHO, NULL, NULL )
    )
  ); 

  /*======================================================================
  // Get right hand side vector.
  //======================================================================*/

  p_rhs = Libnucnet__Zone__computeFlowVector( p_zone );

  p_abund =
    Libnucnet__Zone__getAbundances( p_zone );

  v_old =
    gsl_vector_subvector(
      p_user_data->pYOld,
      i_zone * WnMatrix__get_gsl_vector_size( p_abund ),
      WnMatrix__get_gsl_vector_size( p_abund )
    );

  gsl_vector_sub(
    p_abund, &v_old.vector
  );

  gsl_vector_scale( p_abund, 1. / p_user_data->dDt );

  gsl_vector_sub( p_rhs, p_abund );

  /*======================================================================
  // Insert into large vector.
  //======================================================================*/

  v_rhs =
    gsl_vector_subvector(
      p_user_data->pRhs,
      i_zone * WnMatrix__get_gsl_vector_size( p_rhs ),
      WnMatrix__get_gsl_vector_size( p_rhs )
    );

  gsl_vector_memcpy( &v_rhs.vector, p_rhs );

  gsl_vector_free( p_rhs );

  /*========================================================================
  // Get the Jacobian Matrix.
  //======================================================================*/

  p_matrix = Libnucnet__Zone__computeJacobianMatrix( p_zone );

  /*========================================================================
  // Insert single-zone matrix into full matrix.
  //======================================================================*/

  WnMatrix__insertMatrix(
    p_user_data->pMatrix,
    p_matrix,
    i_zone * WnMatrix__get_gsl_vector_size( p_abund ) + 1,
    i_zone * WnMatrix__get_gsl_vector_size( p_abund ) + 1
  );

  /*========================================================================
  // Free memory allocated for single-zone matrix and abundances. 
  //======================================================================*/

  gsl_vector_free( p_abund );
  WnMatrix__free( p_matrix );

  return 1;

}

/*##############################################################################
// set_mixing_terms().
//############################################################################*/

int
set_mixing_terms( Libnucnet__Zone *p_zone, user_data *p_user_data )
{

  Libnucnet__Zone *p_previous, *p_next;
  gsl_vector *p_abundances;
  gsl_vector_view v_view;
  size_t i, i_species, i_zone;
  char s_label[256];

  i_species =
    Libnucnet__Nuc__getNumberOfSpecies(
      Libnucnet__Net__getNuc(
        Libnucnet__Zone__getNet( p_zone )
      )
    );

  i_zone = (size_t) atol( Libnucnet__Zone__getLabel( p_zone, 1 ) );

  /*============================================================================
  // Add mixing to/from previous zone.
  //==========================================================================*/

  if( i_zone != 0 ) {

    sprintf( s_label, "%lu", (unsigned long) i_zone - 1 );

    v_view =
      gsl_vector_subvector(
        p_user_data->pRhs,
        i_zone * i_species,
        i_species
      );
  
    p_previous =
      Libnucnet__getZoneByLabels(
        p_user_data->pLibnucnet,
        s_label,
        "0",
        "0"
      );

    for( i = 0; i < i_species; i++ )
      WnMatrix__assignElement(
        p_user_data->pMatrix,
        i_zone * i_species + i + 1,
        ( i_zone - 1 ) * i_species + i + 1,
        -1. / p_user_data->dMixTime
      );

    p_abundances = Libnucnet__Zone__getAbundances( p_previous );

    gsl_vector_scale( p_abundances, 1. / p_user_data->dMixTime );

    gsl_vector_add( &v_view.vector, p_abundances );

    gsl_vector_free( p_abundances );

    for( i = 0; i < i_species; i++ )
      WnMatrix__assignElement(
        p_user_data->pMatrix,
        i_zone * i_species + i + 1,
        i_zone * i_species + i + 1,
        1. / p_user_data->dMixTime
      );

    p_abundances = Libnucnet__Zone__getAbundances( p_zone );

    gsl_vector_scale( p_abundances, -1. / p_user_data->dMixTime );

    gsl_vector_add( &v_view.vector, p_abundances );

    gsl_vector_free( p_abundances );

  }

  /*============================================================================
  // Add mixing to/from next zone.
  //==========================================================================*/

  if( i_zone < Libnucnet__getNumberOfZones( p_user_data->pLibnucnet ) - 1 ) {

    sprintf( s_label, "%lu", (unsigned long) i_zone + 1 );

    v_view =
      gsl_vector_subvector(
        p_user_data->pRhs,
        i_zone * i_species,
        i_species
      );
  
    p_next =
      Libnucnet__getZoneByLabels(
        p_user_data->pLibnucnet,
        s_label,
        "0",
        "0"
      );

    for( i = 0; i < i_species; i++ )
      WnMatrix__assignElement(
        p_user_data->pMatrix,
        i_zone * i_species + i + 1,
        ( i_zone + 1 ) * i_species + i + 1,
        -1. / p_user_data->dMixTime
      );

    p_abundances = Libnucnet__Zone__getAbundances( p_next );

    gsl_vector_scale( p_abundances, 1. / p_user_data->dMixTime );

    gsl_vector_add( &v_view.vector, p_abundances );

    gsl_vector_free( p_abundances );

    for( i = 0; i < i_species; i++ )
      WnMatrix__assignElement(
        p_user_data->pMatrix,
        i_zone * i_species + i + 1,
        i_zone * i_species + i + 1,
        1. / p_user_data->dMixTime
      );

    p_abundances = Libnucnet__Zone__getAbundances( p_zone );

    gsl_vector_scale( p_abundances, -1. / p_user_data->dMixTime );

    gsl_vector_add( &v_view.vector, p_abundances );

    gsl_vector_free( p_abundances );

  }

  return 1;

}

/*##############################################################################
// update_abundances().
//############################################################################*/

int
update_abundances( Libnucnet__Zone *p_zone, user_data *p_user_data )
{

  size_t i, i_zone, i_species;
  double d_y, d_checkT;
  gsl_vector *p_new_abundances;
  gsl_vector_view v_view;

  i_zone = (size_t) atol( Libnucnet__Zone__getLabel( p_zone, 1 ) );

  i_species =
    Libnucnet__Nuc__getNumberOfSpecies(
      Libnucnet__Net__getNuc(
        Libnucnet__Zone__getNet( p_zone )
      )
    );

  p_new_abundances = Libnucnet__Zone__getAbundances( p_zone );

  v_view =
    gsl_vector_subvector(
       p_user_data->pSolution,
       i_species * i_zone,
       i_species
    );

  gsl_vector_add( p_new_abundances, &v_view.vector );

  Libnucnet__Zone__updateAbundances( p_zone, p_new_abundances );

  for( i = 0; i < WnMatrix__get_gsl_vector_size( p_new_abundances ); i++ ) {
    d_y = gsl_vector_get( p_new_abundances, i );
    if( d_y > D_Y_MIN ) {
      d_checkT = fabs( gsl_vector_get( &v_view.vector, i ) / d_y );
      if( d_checkT > p_user_data->dCheck )
        p_user_data->dCheck = d_checkT;
    }

  }

  gsl_vector_free( p_new_abundances );

  return 1;

}

/*##############################################################################
// update_abundance_changes().
//############################################################################*/

int
update_abundance_changes( Libnucnet__Zone *p_zone, gsl_vector *p_changes )
{

  size_t i_zone, i_species;
  gsl_vector_view v_view;

  i_zone = (size_t) atol( Libnucnet__Zone__getLabel( p_zone, 1 ) );

  i_species =
    Libnucnet__Nuc__getNumberOfSpecies(
      Libnucnet__Net__getNuc(
        Libnucnet__Zone__getNet( p_zone )
      )
    );

  v_view =
    gsl_vector_subvector(
       p_changes,
       i_zone * i_species,
       i_species
    );

  Libnucnet__Zone__updateAbundanceChanges( p_zone, &v_view.vector );

  return 1;

}

/*##############################################################################
// check_zone_for_t9_and_rho().
//############################################################################*/

int
check_zone_for_t9_and_rho( Libnucnet__Zone *p_zone, int *p_flag )
{

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if(
    !Libnucnet__Zone__getProperty( p_zone, T9, NULL, NULL )
  )
  {
    fprintf(
      stderr,
      "Zone with labels %s %s %s does not have property t9.\n",
      Libnucnet__Zone__getLabel( p_zone, 1 ),
      Libnucnet__Zone__getLabel( p_zone, 2 ),
      Libnucnet__Zone__getLabel( p_zone, 3 )
    );

    *p_flag = 0;
  }

  if(
    !Libnucnet__Zone__getProperty( p_zone, RHO, NULL, NULL )
  )
  {
    fprintf(
      stderr,
      "Zone with labels %s %s %s does not have property rho.\n",
      Libnucnet__Zone__getLabel( p_zone, 1 ),
      Libnucnet__Zone__getLabel( p_zone, 2 ),
      Libnucnet__Zone__getLabel( p_zone, 3 )
    );

    *p_flag = 0;
  }

  return 1;

}

/*##############################################################################
// set_zone_for_mixing_only().
//############################################################################*/

int
set_zone_for_mixing_only( Libnucnet__Zone *p_zone, void *p_data )
{

  if( p_data )
  {
    fprintf( stderr, "Routine should have no extra data.\n" );
    exit( EXIT_FAILURE );
  }

  Libnucnet__Zone__updateProperty( p_zone, T9, NULL, NULL, "0" );

  return 1;

}

/*##############################################################################
// initialize_zone().
//############################################################################*/

int
initialize_zone( Libnucnet__Zone *p_zone, void *p_data )
{

  Libnucnet__Reac *p_my_duplicates;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if( p_data ) {
    fprintf( stderr, "Routine takes no extra data!\n" );
    exit( EXIT_FAILURE );
  }

  /*============================================================================
  // Remove duplicate reactions.
  //==========================================================================*/

  p_my_duplicates =
    Libnucnet__Reac__getDuplicateReactions(
      Libnucnet__Net__getReac(
        Libnucnet__Zone__getNet( p_zone )
      )
    );

  Libnucnet__Reac__iterateReactions(
    p_my_duplicates,
    (Libnucnet__Reaction__iterateFunction) remove_duplicate,
    Libnucnet__Zone__getNet( p_zone )
  );

  Libnucnet__Reac__free( p_my_duplicates );

  return 1;

}

/*##############################################################################
// zone_compare().
//############################################################################*/

int
zone_compare( const Libnucnet__Zone *p_zone1, const Libnucnet__Zone *p_zone2 )
{

  if(
    atol( Libnucnet__Zone__getLabel( p_zone1, 1 ) ) <
    atol( Libnucnet__Zone__getLabel( p_zone2, 1 ) )
  )
    return -1;
  else if(
    atol( Libnucnet__Zone__getLabel( p_zone1, 1 ) ) >
    atol( Libnucnet__Zone__getLabel( p_zone2, 1 ) )
  )
    return 1;
  else
    return 0;

}

/*##############################################################################
// add_copy_of_net_view().
//############################################################################*/

int
add_copy_of_net_view(
  Libnucnet__Zone * p_zone,
  Libnucnet__NetView * p_view
)
{

  return
    Libnucnet__Zone__updateNetView(
      p_zone,
      "",
      "",
      NULL,
      Libnucnet__NetView__copy( p_view )
    ); 

}
