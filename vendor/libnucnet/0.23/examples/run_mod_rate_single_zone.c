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
//       a new Libnucnet structure from an input xml file, evolve the
//       abundances for the input temperature, density, and expansion timescale
//       over the input duration but with the forward and reverse rates
//       for a particular reaction multiplied by a constant factor at all
//       temperatures, and free the allocated memory.
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

#define D_DT0          1.e-15  /* Initial timestep */
#define D_MIN          1.e-4   /* Newton-Raphson convergence criterion */
#define D_REG_T        0.15    /* Largest fraction change in dt */
#define D_REG_Y        0.15    /* Largest change in y for timestep */
#define D_Y_MIN        1.e-10  /* Smallest y for convergence */
#define D_Y_MIN_PRINT  1.e-30  /* Smallest y to print out */
#define D_Y_MIN_DT     1.e-10  /* Smallest y for dt update */

#define I_ITMAX        10      /* Maximum number of Newton-Raphson iterations */
#define I_SOLVER       ARROW   /* Solver type: ARROW or GSL */

#define T9             "t9"    /* String for temperature property */
#define RHO            "rho"   /* String for density property */

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

int
evolve_system(
  Libnucnet__Zone *, double, const char *, double
);
void print_abunds( Libnucnet__Zone *, double, double );
int print_abundance( Libnucnet__Species *, Libnucnet__Zone * );
int
sort_function(
  const Libnucnet__Species *,
  const Libnucnet__Species *
);

/*##############################################################################
// main().
//############################################################################*/

int main( int argc, char * argv[] ) {

  int i_step;
  double d_tau, d_tend, d_t, d_dt;
  Libnucnet *p_my_nucnet;
  Libnucnet__Zone *p_zone;
  Libnucnet__Reac *p_my_duplicates;
  char s_property[32];

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if ( argc < 9 || argc > 11 ) {
    fprintf(
      stderr,
      "\nUsage: %s in_file t9 rho tau t_end m_step reaction factor xpath_nuc xpath_reac\n\n",
      argv[0]
    );
    fprintf(
      stderr, "  in_file = input single zone data xml filename\n\n"
    );
    fprintf(
      stderr, "  t9 = input temperature\n\n"
    );
    fprintf(
      stderr, "  rho = input density (in g/cc)\n\n"
    );
    fprintf(
      stderr, "  tau = density e-folding time scale (<= 0 for static)\n\n"
    );
    fprintf(
      stderr, "  t_end = duration to evolve system\n\n"
    );
    fprintf(
      stderr, "  m_step = frequency of steps to print out\n\n"
    );
    fprintf(
      stderr, "  reaction = string giving reaction to vary\n\n"
    );
    fprintf(
      stderr, "  factor = factor by which to vary reaction\n\n"
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
    if( !Libnucnet__is_valid_input_xml( argv[1] ) ) {
      fprintf( stderr, "Not valid libnucnet input!\n" );
      return EXIT_FAILURE;
    }
  }

  /*============================================================================
  // Read and store input.
  //==========================================================================*/

  if( argc == 9 ) {
    p_my_nucnet = Libnucnet__new_from_xml( argv[1], NULL, NULL, NULL );
  } else if( argc == 10 ) {
    p_my_nucnet = Libnucnet__new_from_xml( argv[1], argv[9], NULL, NULL );
  } else {
    p_my_nucnet = Libnucnet__new_from_xml( argv[1], argv[9], argv[10], NULL );
  }

  /*============================================================================
  // Register the user-supplied rate functions.
  //==========================================================================*/

  register_my_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

  /*============================================================================
  // Initialize time.
  //==========================================================================*/

  d_dt = D_DT0;
  d_t = 0.0;

  /*============================================================================
  // Get zone.
  //==========================================================================*/

  p_zone =
    Libnucnet__getZoneByLabels( p_my_nucnet, "0", "0", "0" );

  if( !p_zone ) {
    fprintf( stderr, "Single zone not found!\n" );
    return EXIT_FAILURE;
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
  // Set temperature and density.
  //==========================================================================*/

  Libnucnet__Zone__updateProperty( p_zone, T9, NULL, NULL, argv[2] );
  
  Libnucnet__Zone__updateProperty( p_zone, RHO, NULL, NULL, argv[3] );

  /*============================================================================
  // Set timescale and end time.
  //==========================================================================*/

  d_tau = atof( argv[4] );
  d_tend = atof( argv[5] );

  /*============================================================================
  // Evolve network while t < final t.
  //==========================================================================*/

  i_step = 1;

  while ( d_t < d_tend ) {

    d_t += d_dt;

  /*============================================================================
  // Get temperature and density.
  //==========================================================================*/

    if( d_tau > 0. ) {
      sprintf(
        s_property, 
        "%g",
        atof( argv[3] ) * exp( -d_t / d_tau )
      );
      Libnucnet__Zone__updateProperty( p_zone, RHO, NULL, NULL, s_property );
      sprintf(
        s_property, 
        "%g",
        atof( argv[2] ) * exp( -d_t / ( 3. * d_tau ) )
      );
      Libnucnet__Zone__updateProperty( p_zone, T9, NULL, NULL, s_property );
    }

  /*============================================================================
  // Evolve abundances.
  //==========================================================================*/

    if( evolve_system( p_zone, d_dt, argv[7], atof( argv[8] ) )
        != 0
    )
      return EXIT_FAILURE;

  /*============================================================================
  // Print out abundances.
  //==========================================================================*/

  if( ( i_step % atoi( argv[6] ) ) == 0 || d_t >= d_tend )
    print_abunds( p_zone, d_t, d_dt );

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

  Libnucnet__free( p_my_nucnet );
  return EXIT_SUCCESS;

}

/*##############################################################################
// evolve_system()
//############################################################################*/

int
evolve_system( 
  Libnucnet__Zone *p_zone,
  double d_dt,
  const char *s_reaction,
  double d_factor 
) {

  Libnucnet__Reaction *p_reaction;
  WnMatrix *p_matrix; 
  WnMatrix__Arrow *p_arrow;
  size_t i, i_species, i_iter;
  gsl_vector *p_y_old, *p_rhs, *p_sol, *p_work;
  double d_forward, d_reverse, d_check, d_checkT;
  double d_rhoe;
  
  /*============================================================================
  // User-supplied data structures.  See screen.h and coul_corr.h for
  // definitions.
  //==========================================================================*/

  user_screening_data my_screening_data;
  user_coul_corr_data my_corr_data = {-0.9052, 0.6322};

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
  // Get reaction to be modified.
  //==========================================================================*/

  p_reaction =
    Libnucnet__Reac__getReactionByString(
      Libnucnet__Net__getReac(
        Libnucnet__Zone__getNet( p_zone )
      ),
      s_reaction
    );

  if( !p_reaction ) {
    fprintf( stderr, "No such reaction!\n" );
    return 1;
  }

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
    // Compute rates.
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
    // Modify reaction.
    //------------------------------------------------------------------------*/

    Libnucnet__Zone__getRatesForReaction(
      p_zone, p_reaction, &d_forward, &d_reverse
    );

    d_forward *= d_factor;
    d_reverse *= d_factor;

    Libnucnet__Zone__updateRatesForReaction(
      p_zone, p_reaction, d_forward, d_reverse
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
     printf( "%5u%5u%14.4e%14.4e\n",
       Libnucnet__Species__getZ( p_species ),
       Libnucnet__Species__getA( p_species ),
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
