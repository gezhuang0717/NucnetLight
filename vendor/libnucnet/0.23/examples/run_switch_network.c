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
//       over the input duration but switch to use a network only including
//       weak decays and neutron captures below a certain temperature,
//       and finally free the allocated memory.
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

#define D_DT0          1.e-10  /* Initial time step */
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

#define VALIDATE       "no"

/*##############################################################################
// Screening.  "no" = no screening, "yes" = screening.
//############################################################################*/

#define USE_SCREENING  "yes"

/*##############################################################################
// Restricted network reaction XPath.  The restriction that the weak decay
// not include h1 as a reactant prevents the network from using the initial
// PP reaction h1 + h1 -> h2 + positron + neutrino_e.
//############################################################################*/

#define S_REAC_NET_XPATH \
   "[ \
      (product = 'electron' or product = 'positron' and not(reactant = 'h1')) \
      or \
      (reactant = 'n' and product = 'gamma') \
   ]"

#define D_T9_NET_RESTRICT  0.1

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

/*##############################################################################
// main().
//############################################################################*/

int main( int argc, char * argv[] ) {

  int i_step;
  double d_tau, d_tend, d_t, d_dt;
  Libnucnet *p_my_nucnet = NULL;
  Libnucnet__NetView * p_view;
  Libnucnet__Zone *p_zone = NULL;
  Libnucnet__Reac *p_duplicates;
  char s_property[32];


  /*============================================================================
  // Check input.
  //==========================================================================*/

  if ( argc < 10 || argc > 12 ) {
    fprintf(
      stderr,
      "\nUsage: %s in_file label_1 label_2 label_3 t9 rho tau t_end m_step xpath_nuc xpath_reac\n\n",
      argv[0]
    );
    fprintf(
      stderr, "  in_file = input single zone data xml filename\n\n"
    );
    fprintf(
      stderr, "  label_1 = first label of zone\n\n"
    );
    fprintf(
      stderr, "  label_2 = second label of zone\n\n"
    );
    fprintf(
      stderr, "  label_3 = third label of zone\n\n"
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

  if( argc == 10 ) {
    p_my_nucnet = Libnucnet__new_from_xml( argv[1], NULL, NULL, NULL );
  } else if( argc == 11 ) {
    p_my_nucnet = Libnucnet__new_from_xml( argv[1], argv[10], NULL, NULL );
  } else {
    p_my_nucnet = Libnucnet__new_from_xml( argv[1], argv[10], argv[11], NULL );
  }

  d_tau = atof( argv[7] );
  d_tend = atof( argv[8] );

  /*============================================================================
  // Register user-supplied rate functions.
  //==========================================================================*/

  register_my_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

  /*============================================================================
  // Since we will allocate data for some of the user rate functions,
  // we set the rate function data deallocators.
  //==========================================================================*/

  set_user_data_deallocators(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

  /*============================================================================
  // Initialize time.
  //==========================================================================*/

  d_dt = D_DT0;
  d_t = 0.0;

  /*============================================================================
  // Get the zone.
  //==========================================================================*/

  p_zone =
    Libnucnet__getZoneByLabels( p_my_nucnet, argv[2], argv[3], argv[4] );

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
  // Set temperature and density.
  //==========================================================================*/

  Libnucnet__Zone__updateProperty( p_zone, T9, NULL, NULL, argv[5] );
  
  Libnucnet__Zone__updateProperty( p_zone, RHO, NULL, NULL, argv[6] );

  /*============================================================================
  // Evolve network while t < final t.
  //==========================================================================*/

  i_step = 0;

  while ( d_t < d_tend ) {

    d_t += d_dt;

  /*============================================================================
  // Get temperature and density.
  //==========================================================================*/

    if( d_tau > 0. )
    { 
      sprintf(
        s_property, 
        "%g",
        atof( argv[6] ) * exp( -d_t / d_tau )
      );
      Libnucnet__Zone__updateProperty( p_zone, RHO, NULL, NULL, s_property );
      sprintf(
        s_property, 
        "%g",
        atof( argv[5] ) * exp( -d_t / ( 3. * d_tau ) )
      );
      Libnucnet__Zone__updateProperty( p_zone, T9, NULL, NULL, s_property );
    }

  /*============================================================================
  // Switch to use only weak decays and (n,g) and (g,n) below T9 = 0.1.
  // We create a view, store it in the zone for future lookup, and store a
  // copy as the evolution network.
  //==========================================================================*/

    if(
      atof( Libnucnet__Zone__getProperty( p_zone, T9, NULL, NULL ) )
      <
      D_T9_NET_RESTRICT
    ) 
    {

      p_view =
        Libnucnet__Zone__getNetView(
          p_zone,
          "",
          S_REAC_NET_XPATH,
          NULL
        );

      if( !p_view || Libnucnet__NetView__wasNetUpdated( p_view ) )
      {

        p_view =
          Libnucnet__NetView__new(
            Libnucnet__Zone__getNet( p_zone ),
            "",
            S_REAC_NET_XPATH
          );

        Libnucnet__Zone__updateNetView(
          p_zone,
          "",
          S_REAC_NET_XPATH,
          NULL,
          p_view
        );

        Libnucnet__Zone__updateNetView(
          p_zone,
          EVOLUTION_NETWORK,
          NULL,
          NULL,
          Libnucnet__NetView__copy( p_view )
        );

      }

    }

  /*============================================================================
  // Evolve abundances.
  //==========================================================================*/

    if ( evolve_system( p_zone, d_dt ) != 0 ) {

      return EXIT_FAILURE;

    }

  /*============================================================================
  // Print out abundances.
  //==========================================================================*/

  if( ( i_step % atoi( argv[9] ) ) == 0 || d_t >= d_tend )
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
  // Clear net views and rate values stored in p_zone.  It is not necessary
  // to do this since they are cleared when Libnucnet__free clears the zone.
  // We do it for illustration.
  //==========================================================================*/

  if( !Libnucnet__Zone__clearNetViews( p_zone ) )
  {
    fprintf( stderr, "Couldn't clear net views.\n" );
  }

  if( !Libnucnet__Zone__clearRates( p_zone ) )
  {
    fprintf( stderr, "Couldn't clear rate values.\n" );
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
  Libnucnet__Zone *p_zone, double d_dt 
) {

  WnMatrix *p_matrix; 
  WnMatrix__Arrow *p_arrow;
  size_t i, i_species, i_iter;
  gsl_vector *p_y_old, *p_rhs, *p_sol, *p_work;
  double *p_rhoe, d_check, d_checkT;
  user_screening_data my_screening_data;
  user_coul_corr_data my_corr_data = {-0.9052,0.6322};
  
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

    p_rhoe = ( double * ) malloc( sizeof( double ) );

    *p_rhoe =
      atof(
        Libnucnet__Zone__getProperty( p_zone, RHO, NULL, NULL )
      ) *
      Libnucnet__Zone__computeZMoment( p_zone, 1 );
  
    Libnucnet__Zone__updateDataForUserRateFunction(
      p_zone,
      CF88_WEAK_FIT,
      p_rhoe
    );

    /*--------------------------------------------------------------------------
    // Set screening and coulomb correction functions. Update screening data.
    // See screen.h and coul_corr.h for definitions.
    //------------------------------------------------------------------------*/

    if( strcmp( USE_SCREENING, "yes" ) == 0 )
    {
  
      my_screening_data.dYe2 =
        Libnucnet__Zone__computeZMoment( p_zone, 2 );

      Libnucnet__Zone__setScreeningFunction(
        p_zone,
        (Libnucnet__Net__screening_function) my_screening_function,
        &my_screening_data
      );
    
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

    p_matrix = Libnucnet__Zone__computeJacobian( p_zone );
      WnMatrix__scaleMatrix( p_matrix, -1. );
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
   
