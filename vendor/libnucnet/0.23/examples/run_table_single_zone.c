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
//       a new Libnucnet structure from input xml and ascii files, evolve the
//       abundances for the input temperature, density, and expansion timescale
//       over the input duration, and free the allocated memory.
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
// Enumeration.
//############################################################################*/

enum solvers { ARROW, GSL };

/*##############################################################################
// Define some parameters.
//############################################################################*/

#define BUF            32      /* String buffer size */
#define D_DT0          1.e-10  /* Initial time step */
#define D_EPS          0.05    /* The maximum change in t9 or rho */
#define D_MIN          1.e-4   /* Newton-Raphson convergence criterion */
#define D_REG_T        0.15    /* Time step change regulator for dt update */
#define D_REG_Y        0.15    /* Abundance change regulator for dt update */
#define D_Y_MIN        1.e-10  /* Smallest y for convergence */
#define D_Y_MIN_DT     1.e-10  /* Smallest y for dt update */
#define D_Y_MIN_PRINT  1.e-30  /* Smallest y to print out */
#define I_ITMAX        10      /* Maximum number of Newton-Raphson iterations */
#define I_SOLVER       ARROW   /* Solver type: ARROW or GSL */

/*##############################################################################
// Define some strings.
//############################################################################*/

#define T9 "t9"                /* A string for the temperature */
#define RHO "rho"              /* A string for the density */

/*##############################################################################
// Validation.  "no" = no validation, "yes" = validation.
//############################################################################*/

#define VALIDATE       "yes"

/*##############################################################################
// Screening.  "no" = no screening, "yes" = screening.
//############################################################################*/

#define USE_SCREENING  "yes"

/*##############################################################################
// Prototypes.
//############################################################################*/

int evolve_system( Libnucnet__Zone *, double );
void print_abunds( Libnucnet__Zone *, double, double );
int print_abundance( Libnucnet__Species *, Libnucnet__Zone * );
int
sort_function(
  const Libnucnet__Species *,
  const Libnucnet__Species *
);
double interpolate_array( double, gsl_vector *, gsl_vector * );

int
update_temperature_and_density(
  Libnucnet__Zone *,
  gsl_vector *,
  gsl_vector *,
  gsl_vector *,
  double *,
  double *
);

/*##############################################################################
// main().
//############################################################################*/

int main( int argc, char * argv[] ) {

  int i_step;
  size_t i, i_size;
  gsl_vector *p_t, *p_t9, *p_log10_rho;
  double d_tend, d_t, d_dt;
  double d_x1, d_x2, d_x3;
  Libnucnet *p_my_nucnet;
  Libnucnet__Zone *p_zone;
  Libnucnet__Reac *p_duplicates;
  Libnucnet__Species *p_species;
  char s_species[BUF];
  FILE *p_file;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if ( argc < 6 || argc > 8 ) {
    fprintf(
      stderr,
      "\nUsage: %s net_file zone_ascii zone_file t_end m_step xpath_nuc xpath_reac\n\n",
      argv[0]
    );
    fprintf(
      stderr, "  net_file = nuclear network data xml filename\n\n"
    );
    fprintf(
      stderr, "  zone_ascii = input single zone ascii filename\n\n"
    );
    fprintf(
      stderr, "  zone_file = input single zone mass fraction filename\n\n"
    );
    fprintf(
      stderr, "  t_end = duration to evolve system\n\n"
    );
    fprintf(
      stderr, "  m_step = frequency of steps to print out\n\n"
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

  /*============================================================================
  // Validate input file.
  //==========================================================================*/

  if( strcmp( VALIDATE, "yes" ) == 0 )
  {
    if( !Libnucnet__Net__is_valid_input_xml( argv[1] ) ) {
      fprintf( stderr, "Not valid libnucnet input!\n" );
      return EXIT_FAILURE;
    }
  }

  /*============================================================================
  // Read and store input.
  //==========================================================================*/

  p_my_nucnet = Libnucnet__new();

  if( argc == 6 )
    Libnucnet__Net__updateFromXml(
      Libnucnet__getNet( p_my_nucnet ),
      argv[1],
      NULL,
      NULL
    );  
  else if( argc == 7 )
    Libnucnet__Net__updateFromXml(
      Libnucnet__getNet( p_my_nucnet ),
      argv[1],
      argv[6],
      NULL
    );  
  else
    Libnucnet__Net__updateFromXml(
      Libnucnet__getNet( p_my_nucnet ),
      argv[1],
      argv[6],
      argv[7]
    );  

  /*============================================================================
  // Register user-supplied rate functions.
  //==========================================================================*/

  register_my_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

  /*============================================================================
  // Create zone and add to nuclear network.
  //==========================================================================*/

  p_zone =
    Libnucnet__Zone__new( Libnucnet__getNet( p_my_nucnet ), "0","0", "0" );

  Libnucnet__addZone( p_my_nucnet, p_zone );

  /*============================================================================
  // Open file thermodynamics file.
  //==========================================================================*/

  p_file = fopen( argv[2], "r" );

  if( !p_file )
  {
    fprintf( stderr, "Couldn't open file.\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Count rows.
  //==========================================================================*/

  i_size = 0;
  do
  {
    fscanf( p_file, "%lf  %lf  %lf\n", &d_x1, &d_x2, &d_x3 );
    i_size++;
  } while( !feof( p_file ) );
  rewind( p_file );

  /*============================================================================
  // Assign arrays.
  //==========================================================================*/

  p_t = gsl_vector_alloc( i_size );
  p_t9 = gsl_vector_alloc( i_size );
  p_log10_rho = gsl_vector_alloc( i_size );

  for( i = 0; i < i_size; i++ )
  {
    fscanf( p_file, "%lf  %lf  %lf\n", &d_x1, &d_x2, &d_x3 );
    gsl_vector_set( p_t, i, d_x1 );
    gsl_vector_set( p_t9, i, d_x2 );
    gsl_vector_set( p_log10_rho, i, log10( d_x3 ) );
  }

  fclose( p_file );

  /*============================================================================
  // Open file mass fractions file.
  //==========================================================================*/

  p_file = fopen( argv[3], "r" );

  if( !p_file )
  {
    fprintf( stderr, "Couldn't open file.\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Read in mass fraction data.
  //==========================================================================*/

  do
  {
    fscanf( p_file, "%s  %lf\n", s_species, &d_x1 );

    p_species =
      Libnucnet__Nuc__getSpeciesByName(
        Libnucnet__Net__getNuc(
          Libnucnet__getNet( p_my_nucnet )
        ),
        s_species
      );

    Libnucnet__Zone__updateSpeciesAbundance(
      p_zone,
      p_species,
      d_x1 / Libnucnet__Species__getA( p_species )
    );
      
  } while( !feof( p_file ) );

  fclose( p_file );

  /*============================================================================
  // Initialize time parameters.
  //==========================================================================*/

  d_dt = D_DT0;
  d_t = 0.0;
  d_tend = atof( argv[4] );

  /*============================================================================
  // Get the zone.
  //==========================================================================*/

  p_zone =
    Libnucnet__getZoneByLabels( p_my_nucnet, "0", "0", "0" );

  /*============================================================================
  // Remove duplicate reactions.
  //==========================================================================*/

  p_duplicates =
    Libnucnet__Reac__getDuplicateReactions(
      Libnucnet__Net__getReac(
        Libnucnet__Zone__getNet( p_zone )
      )
    );

  Libnucnet__Reac__iterateReactions(
    p_duplicates,
    (Libnucnet__Reaction__iterateFunction) remove_duplicate,
    Libnucnet__Zone__getNet( p_zone )
  );

  Libnucnet__Reac__free( p_duplicates );

  /*============================================================================
  // Sort the nuclei if using the arrow solver.
  //==========================================================================*/

  if( I_SOLVER == ARROW )
  {

    Libnucnet__Nuc__setSpeciesCompareFunction(
      Libnucnet__Net__getNuc( Libnucnet__Zone__getNet( p_zone ) ),
      (Libnucnet__Species__compare_function) sort_function
    );

    Libnucnet__Nuc__sortSpecies(
      Libnucnet__Net__getNuc( Libnucnet__Zone__getNet( p_zone ) ) 
    );

  }

  /*============================================================================
  // Evolve network while t < final t.
  //==========================================================================*/

  i_step = 0;

  while ( d_t < d_tend ) {

    d_t += d_dt;

  /*============================================================================
  // Get temperature and density.
  //==========================================================================*/

    update_temperature_and_density(
      p_zone,
      p_t,
      p_t9,
      p_log10_rho,
      &d_t,
      &d_dt
    );

  /*============================================================================
  // Evolve abundances.
  //==========================================================================*/

    if ( evolve_system( p_zone, d_dt ) != 0 ) {

      return EXIT_FAILURE;

    }

  /*============================================================================
  // Print out abundances.
  //==========================================================================*/

  if( ( i_step % atoi( argv[5] ) ) == 0 || d_t >= d_tend ) {
    print_abunds( p_zone, d_t, d_dt );
  }

  /*============================================================================
  // Update timestep.
  //==========================================================================*/

    Libnucnet__Zone__updateTimeStep(
      p_zone,
      &d_dt,
      D_REG_T,
      D_REG_Y,
      D_Y_MIN_DT
    );

    if ( d_t + d_dt > d_tend ) {

      d_dt = d_tend - d_t;

    }

    i_step++;

  }  

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  gsl_vector_free( p_t );
  gsl_vector_free( p_t9 );
  gsl_vector_free( p_log10_rho );

  Libnucnet__free( p_my_nucnet );
  return EXIT_SUCCESS;

}

/*##############################################################################
// evolve_system()
//############################################################################*/

int
evolve_system( 
  Libnucnet__Zone *p_zone, double d_dt 
) {

  WnMatrix *p_matrix; 
  WnMatrix__Arrow *p_arrow;
  size_t i, i_species, i_iter;
  gsl_vector *p_y_old, *p_rhs, *p_sol, *p_work;
  double d_check, d_checkT;
  double d_rhoe;
  
  /*============================================================================
  // User-supplied data structures.  See screen.h and coul_corr.h for
  // definitions.
  //==========================================================================*/

  user_screening_data my_screening_data;
  user_coul_corr_data my_corr_data = {-0.9052,0.6322};

  /*============================================================================
  // Get number of species.
  //==========================================================================*/

  i_species =
    Libnucnet__Nuc__getNumberOfSpecies(
      Libnucnet__Net__getNuc(
        Libnucnet__Zone__getNet( p_zone )
      )
    );

  /*============================================================================
  // Save the old abundances.
  //==========================================================================*/

  p_y_old = Libnucnet__Zone__getAbundances( p_zone );

  /*============================================================================
  // Newton-Raphson Iterations.
  //==========================================================================*/

  for( i_iter = 0; i_iter < I_ITMAX; i_iter++ ) {

    /*--------------------------------------------------------------------------
    // Update data for user-rate functions.
    //------------------------------------------------------------------------*/

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
      
    /*--------------------------------------------------------------------------
    // Set screening, if desired.
    //------------------------------------------------------------------------*/

    if( strcmp( USE_SCREENING, "yes" ) == 0 ) {

      /*........................................................................
      // Set screening data and function.
      //......................................................................*/

      my_screening_data.dYe2 =
        Libnucnet__Zone__computeZMoment( p_zone, 2 );

      Libnucnet__Zone__setScreeningFunction(
        p_zone,
        (Libnucnet__Net__screening_function) my_screening_function,
        &my_screening_data
      );

      /*........................................................................
      // Set correction factor function.
      //......................................................................*/

      Libnucnet__Zone__setNseCorrectionFactorFunction(
        p_zone,
        (Libnucnet__Species__nseCorrectionFactorFunction) my_coulomb_correction,
        &my_corr_data
      );

    }

    /*--------------------------------------------------------------------------
    // Compute rates.  For iterations after the first, only update rates
    // for those reactions for which the rate function has been tagged for
    // updating.
    //------------------------------------------------------------------------*/

    Libnucnet__Zone__computeRates(
      p_zone,
      atof(
        Libnucnet__Zone__getProperty( p_zone, T9, NULL, NULL )
      ),
      atof(
        Libnucnet__Zone__getProperty( p_zone, RHO, NULL, NULL )
      )
    ); 

    /*--------------------------------------------------------------------------
    // Get right hand side vector.
    //------------------------------------------------------------------------*/

    p_rhs = Libnucnet__Zone__computeFlowVector( p_zone );

    p_work = Libnucnet__Zone__getAbundances( p_zone );
    gsl_vector_sub( p_work, p_y_old );
    gsl_vector_scale( p_work, 1. / d_dt );
    gsl_vector_sub( p_rhs, p_work );

    gsl_vector_free( p_work );

    /*--------------------------------------------------------------------------
    // Get Jacobian matrix.
    //------------------------------------------------------------------------*/

    p_matrix = Libnucnet__Zone__computeJacobianMatrix( p_zone );
    WnMatrix__addValueToDiagonals( p_matrix, 1.0 / d_dt );

    /*--------------------------------------------------------------------------
    // Solve matrix equation.
    //------------------------------------------------------------------------*/

    if( I_SOLVER == ARROW )
    {
      p_arrow = WnMatrix__getArrow( p_matrix, 3L );
      p_sol = WnMatrix__Arrow__solve( p_arrow, p_rhs );
      WnMatrix__Arrow__free( p_arrow );
    }
    else
      p_sol = WnMatrix__solve( p_matrix, p_rhs );
  
    /*--------------------------------------------------------------------------
    // Update abundances and check for convergence.
    //------------------------------------------------------------------------*/

    d_check = 0.;

    p_work = Libnucnet__Zone__getAbundances( p_zone );

    for( i = 0; i < i_species; i++ ) {
      if( gsl_vector_get( p_work, i )  > D_Y_MIN ) {
        d_checkT =
          fabs( gsl_vector_get( p_sol, i ) / gsl_vector_get( p_work, i ));
        if( d_checkT > d_check ) d_check = d_checkT;
      }
    }
    
    gsl_vector_add( p_work, p_sol );

    Libnucnet__Zone__updateAbundances( p_zone, p_work );

    gsl_vector_free( p_work );

    /*--------------------------------------------------------------------------
    // Free matrix, p_rhs, and p_sol. Remember
    // Libnucnet__Zone__computeJacobianMatrix returns a new matrix and
    // Libnucnet__computeFlowVector and WnMatrix__solve return new gsl_vectors
    // each time they are called.
    //------------------------------------------------------------------------*/

    WnMatrix__free( p_matrix ); 
    gsl_vector_free( p_rhs );
    gsl_vector_free( p_sol );

    /*--------------------------------------------------------------------------
    // Exit iterations if converged.
    //------------------------------------------------------------------------*/

    if( d_check < D_MIN ) break;

  }

  /*==========================================================================
  // Update abundance changes.
  //========================================================================*/

  p_work = Libnucnet__Zone__getAbundances( p_zone );

  gsl_vector_sub( p_work, p_y_old );

  Libnucnet__Zone__updateAbundanceChanges( p_zone, p_work );

  gsl_vector_free( p_work );
  
  /*==========================================================================
  // Free allocated memory and return.
  //========================================================================*/

  gsl_vector_free( p_y_old );

  return 0;

}

/*##############################################################################
// print_abunds()
//############################################################################*/

void
print_abunds(
  Libnucnet__Zone *p_zone, double d_t, double d_dt
) {

  printf(
    "t = %10.4e, dt = %10.4e, t9 = %10.4e, rho (g/cc) = %10.4e\n\n",
    d_t,
    d_dt,
    atof(
      Libnucnet__Zone__getProperty( p_zone, T9, NULL, NULL )
    ),
    atof(
      Libnucnet__Zone__getProperty( p_zone, RHO, NULL, NULL )
    )
  );

  Libnucnet__Nuc__iterateSpecies(
    Libnucnet__Net__getNuc(
      Libnucnet__Zone__getNet( p_zone )
    ),
    (Libnucnet__Species__iterateFunction) print_abundance,
    p_zone
  );

  fprintf(
    stdout,
    "1 - xsum = %e\n\n",
    1. - Libnucnet__Zone__computeAMoment( p_zone, 1 )
  );

}

/*##############################################################################
// print_abundance()
//############################################################################*/

int
print_abundance( Libnucnet__Species *p_species, Libnucnet__Zone *p_zone )
{

  double d_abund;

  if( 
      ( d_abund =
         Libnucnet__Zone__getSpeciesAbundance( p_zone, p_species )
      ) > D_Y_MIN_PRINT
  )
  {
     printf( "%5lu%5lu%14.4e%14.4e\n",
       (unsigned long) Libnucnet__Species__getZ( p_species ),
       (unsigned long) Libnucnet__Species__getA( p_species ),
       d_abund,
       Libnucnet__Zone__getSpeciesAbundanceChange(
         p_zone, p_species
     )
    );
  }

  return 1;

}
/*##############################################################################
// sort_function()
//############################################################################*/

int
sort_function(
  const Libnucnet__Species *p_species1,
  const Libnucnet__Species *p_species2
)
{

  int i;

  if( !strcmp( Libnucnet__Species__getName( p_species1 ), "he4" ) )
    return 1;
  else if( !strcmp( Libnucnet__Species__getName( p_species2 ), "he4" ) )
    return -1;

  if( !strcmp( Libnucnet__Species__getName( p_species1 ), "h1" ) )
    return 1;
  else if( !strcmp( Libnucnet__Species__getName( p_species2 ), "h1" ) )
    return -1;

  if( !strcmp( Libnucnet__Species__getName( p_species1 ), "n" ) )
    return 1;
  else if( !strcmp( Libnucnet__Species__getName( p_species2 ), "n" ) )
    return -1;

  if( 
      Libnucnet__Species__getZ( p_species1 ) <
      Libnucnet__Species__getZ( p_species2 )
  )
    return -1;
  else if(
      Libnucnet__Species__getZ( p_species1 ) >
      Libnucnet__Species__getZ( p_species2 )
  )
    return 1;

  if( 
      Libnucnet__Species__getA( p_species1 ) <
      Libnucnet__Species__getA( p_species2 )
  )
    return -1;
  else if(
      Libnucnet__Species__getA( p_species1 ) >
      Libnucnet__Species__getA( p_species2 )
  )
    return 1;

  i =
    strcmp(
      Libnucnet__Species__getName( p_species1 ),
      Libnucnet__Species__getName( p_species2 )
    );

  if( i == 0 ) {
    return 0;
  } else {
    return GSL_SIGN( i );
  }

}
   
/*##############################################################################
// interpolate_array()
//############################################################################*/

double
interpolate_array(
  double d_t, gsl_vector *p_x, gsl_vector *p_y
)
{

  gsl_interp_accel *acc;
  gsl_spline *spline;
  double d_result;

  if( d_t < WnMatrix__get_gsl_vector_array( p_x )[0] )
    return WnMatrix__get_gsl_vector_array( p_y )[0];

  if(
     d_t >=
       WnMatrix__get_gsl_vector_array(
         p_x
       )[WnMatrix__get_gsl_vector_size( p_x ) - 1]
  )
    return
      WnMatrix__get_gsl_vector_array(
         p_y
       )[WnMatrix__get_gsl_vector_size( p_y ) - 1];

  acc = gsl_interp_accel_alloc();

  spline =
    gsl_spline_alloc(
      gsl_interp_cspline, 
      WnMatrix__get_gsl_vector_size( p_x )
    );

  gsl_spline_init(
    spline,
    WnMatrix__get_gsl_vector_array( p_x ),
    WnMatrix__get_gsl_vector_array( p_y ),
    WnMatrix__get_gsl_vector_size( p_x )
  );

  d_result = gsl_spline_eval( spline, d_t, acc );

  gsl_spline_free( spline );
  gsl_interp_accel_free( acc );

  return d_result;

}

/*##############################################################################
// update_temperature_and_density()
//############################################################################*/

int
update_temperature_and_density(
  Libnucnet__Zone *p_zone,
  gsl_vector *p_t,
  gsl_vector *p_t9,
  gsl_vector *p_log10_rho,
  double *p_time,
  double *p_dt
)
{

  static double d_old_t9 = 0., d_old_rho = 0.;
  static int i_flag = 0;
  double d_t9, d_rho, d_change_t9, d_change_rho;
  char s_property[BUF];
  int i = 0;

  /*==========================================================================
  // Interpolate.  If change in t9 or rho too large, decrease time step by
  // factor of two and interpolate again.
  //========================================================================*/

  do{

     ++i;

     d_t9 = interpolate_array( *p_time, p_t, p_t9 );
     d_rho = pow( 10., interpolate_array( *p_time, p_t, p_log10_rho ) );

     if( i_flag == 0 )
     {

       d_old_t9 = d_t9;
       d_old_rho = d_rho;
       d_change_t9 = 0.;
       d_change_rho = 0.;
       i_flag = 1;

     }
     else
     {

       d_change_t9 = fabs( d_t9 - d_old_t9 ) / d_old_t9;
       d_change_rho = fabs( d_rho - d_old_rho ) / d_old_rho;

     }

     if( d_change_t9 < D_EPS && d_change_rho < D_EPS )
       break;

     *p_time -= *p_dt;
     *p_dt /= 1.1;
     *p_time += *p_dt;
        
  } while( *p_dt > 0. );
  
  sprintf(
    s_property, 
    "%g",
    d_t9
  );
  Libnucnet__Zone__updateProperty( p_zone, T9, NULL, NULL, s_property );

  sprintf(
    s_property, 
    "%g",
    d_rho
  );
  Libnucnet__Zone__updateProperty( p_zone, RHO, NULL, NULL, s_property );

  d_old_t9 = d_t9;
  d_old_rho = d_rho;

  return i;

}
