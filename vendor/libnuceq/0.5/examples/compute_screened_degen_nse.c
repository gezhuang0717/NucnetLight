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
//       Example to demonstrate how to use libnuceq routines to create a
//       an equilibrium structure, compute the default nuclear statistical
//       equilibrium (NSE), NSE with Coulomb correction, and with degenerate,
//       fully relativistic nucleons, and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/
 
#include <Libnuceq.h>
#include "coul_corr.h"

/*##############################################################################
// defines.
//############################################################################*/

#define D_MIN_X_PRINT    1.e-10    /* Minimum mass fraction to print out. */

#define S_NUCEQ_NEUTRON  "n"
#define S_NUCEQ_PROTON   "h1"

/*##############################################################################
// prototypes.
//############################################################################*/

int print_mass_fractions( Libnuceq__Species *, void * );

double
degenerate_function(
  Libnucnet__Species *,
  double,
  double,
  double,
  Libnuceq *
);

double
nse_exp_factor(
  Libnucnet__Nuc *,
  Libnucnet__Species *,
  double,
  double
);

/*##############################################################################
// main().
//############################################################################*/

int main( int argc, char **argv ) {

  typedef struct {
    gsl_vector *pDefault;
    gsl_vector *pCoulombCorrected;
    gsl_vector *pDegenerateNucleons;
  } user_data;

  Libnucnet__Nuc *p_my_nuclei;
  Libnuceq *p_my_equil;
  user_coul_corr_data my_corr_data = {-0.9052,0.6322};
  user_data *p_user_data;
  double d_mun_kT, d_mup_kT;
  double d_mun_kT_corr, d_mup_kT_corr;
  double d_mun_kT_degen, d_mup_kT_degen;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if ( argc < 5 || argc > 6) {
      fprintf(
        stderr, "\nUsage: %s input_file t9 rho Ye nuc_xpath\n", argv[0]
      );
      fprintf(
        stderr, "\n  input_file: name of input xml file\n"
      );
      fprintf(
        stderr, "\n  t9: Temperature in 10^9 K\n"
      );
      fprintf(
        stderr, "\n  rho: mass density in g/cc\n"
      );
      fprintf(
        stderr, "\n  Ye: electron-to-nucleon ratio\n"
      );
      fprintf(
        stderr, "\n  nuc_xpath: XPath expression for nuclides (optional)\n\n"
      );
      return EXIT_FAILURE;
   }

  /*============================================================================
  // Get the nuclide collection.
  //==========================================================================*/

  if( argc == 6 )
    p_my_nuclei = Libnucnet__Nuc__new_from_xml( argv[1], argv[5] );
  else
    p_my_nuclei = Libnucnet__Nuc__new_from_xml( argv[1], NULL );

  /*============================================================================
  // Create equilibrium.
  //==========================================================================*/

  p_my_equil = Libnuceq__new( p_my_nuclei );

  Libnuceq__setYe( p_my_equil, atof( argv[4] ) );

  /*============================================================================
  // Set up work structure.
  //==========================================================================*/
  
  p_user_data = ( user_data * ) malloc( sizeof( user_data ) );

  if( !p_user_data )
  {
    fprintf( stderr, "Couldn't allocate memory.\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Solve the equilibrium for default case.
  //==========================================================================*/

  Libnuceq__computeEquilibrium( p_my_equil, atof( argv[2] ), atof( argv[3] ) );

  p_user_data->pDefault = Libnuceq__getAbundances( p_my_equil );

  d_mun_kT = Libnuceq__getMunkT( p_my_equil );
  d_mup_kT = Libnuceq__getMupkT( p_my_equil );

  /*============================================================================
  // Set correction.
  //==========================================================================*/

  Libnuceq__setNseCorrectionFactorFunction(
    p_my_equil,
    (Libnucnet__Species__nseCorrectionFactorFunction) my_coulomb_correction,
    &my_corr_data
  );

  /*============================================================================
  // Solve the equilibrium for screened case.
  //==========================================================================*/

  Libnuceq__computeEquilibrium( p_my_equil, atof( argv[2] ), atof( argv[3] ) );

  p_user_data->pCoulombCorrected = Libnuceq__getAbundances( p_my_equil );

  d_mun_kT_corr = Libnuceq__getMunkT( p_my_equil );
  d_mup_kT_corr = Libnuceq__getMupkT( p_my_equil );

  /*============================================================================
  // Solve the equilibrium for degenerate, non-screened case.
  //==========================================================================*/
  
  Libnuceq__clearNseCorrectionFactorFunction( p_my_equil );

  Libnuceq__setNseCorrectionFactorFunction(
    p_my_equil,
    (Libnucnet__Species__nseCorrectionFactorFunction) degenerate_function,
    p_my_equil
  );

  Libnuceq__computeEquilibrium( p_my_equil, atof( argv[2] ), atof( argv[3] ) );

  p_user_data->pDegenerateNucleons = Libnuceq__getAbundances( p_my_equil );

  d_mun_kT_degen = Libnuceq__getMunkT( p_my_equil );
  d_mup_kT_degen = Libnuceq__getMupkT( p_my_equil );

  /*============================================================================
  // Print abundances.
  //==========================================================================*/

  fprintf(
    stdout,
    "\nSpecies         X         X (screened)    X (degenerate n,p)\n"
  );  

  fprintf(
    stdout,
    "-------   -------------   -------------   ------------------\n"
  );  

  Libnuceq__iterateSpecies(
    p_my_equil,
    (Libnuceq__Species__iterateFunction) print_mass_fractions,
    p_user_data
  );

  /*============================================================================
  // Print out diagnostics.
  //==========================================================================*/

  fprintf(
    stdout,
    "\nGlobal results for a T9 = %f, rho = %e g/cc, Ye = %f NSE:\n",
    Libnuceq__getT9( p_my_equil ),
    Libnuceq__getRho( p_my_equil ),
    Libnuceq__computeZMoment( p_my_equil, 1 )
  );

  fprintf(
    stdout,
    "\nDefault:\n"
  );

  fprintf(
    stdout,
    "\nNeutron chemical potential / kT = %f\n",
    d_mun_kT
  );

  fprintf(
    stdout,
    "Proton chemical potential / kT = %f\n",
    d_mup_kT
  );

  fprintf(
    stdout,
    "\nScreened:\n"
  );

  fprintf(
    stdout,
    "\nNeutron chemical potential / kT = %f\n",
    d_mun_kT_corr
  );

  fprintf(
    stdout,
    "Proton chemical potential / kT = %f\n",
    d_mup_kT_corr
  );

  fprintf(
    stdout,
    "\nDegenerate:\n"
  );

  fprintf(
    stdout,
    "\nNeutron chemical potential / kT = %f\n",
    d_mun_kT_degen
  );

  fprintf(
    stdout,
    "Proton chemical potential / kT = %f\n",
    d_mup_kT_degen
  );

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  gsl_vector_free( p_user_data->pDefault );
  gsl_vector_free( p_user_data->pCoulombCorrected );
  gsl_vector_free( p_user_data->pDegenerateNucleons );
  free( p_user_data );

  Libnuceq__free( p_my_equil );
  Libnucnet__Nuc__free( p_my_nuclei );

  /*============================================================================
  // Done!
  //==========================================================================*/

  return EXIT_SUCCESS;

}

/*##############################################################################
// print_mass_fractions()
//############################################################################*/

int
print_mass_fractions(
  Libnuceq__Species *p_eq_species,
  void *p_data
)
{

  typedef struct {
    gsl_vector *pDefault;
    gsl_vector *pCoulombCorrected;
    gsl_vector *pDegenerateNucleons;
  } user_data;

  size_t i;
  unsigned int i_a;
  Libnucnet__Species *p_species;
  user_data * p_user_data = ( user_data * ) p_data;
  double d_x, d_x_corr, d_x_degen;

  p_species = Libnuceq__Species__getNucSpecies( p_eq_species );

  i = Libnucnet__Species__getIndex( p_species );

  i_a = Libnucnet__Species__getA( p_species );

  d_x = gsl_vector_get( p_user_data->pDefault, i ) * i_a;
      
  d_x_corr = gsl_vector_get( p_user_data->pCoulombCorrected, i ) * i_a;
      
  d_x_degen = gsl_vector_get( p_user_data->pDegenerateNucleons, i ) * i_a;
      
  if(
     d_x > D_MIN_X_PRINT ||
     d_x_corr > D_MIN_X_PRINT ||
     d_x_degen > D_MIN_X_PRINT
  )
    fprintf(
      stdout,
      "%5s%18.7e%16.7e%16.7e\n",
      Libnucnet__Species__getName( p_species ),
      d_x,
      d_x_corr,
      d_x_degen
    );

  return 1;

}

/*##############################################################################
// degenerate_function().
//############################################################################*/

double
degenerate_function(
  Libnucnet__Species *p_species,
  double d_T9,
  double d_rho,
  double d_ye,
  Libnuceq *p_equil
)
{

  Libstatmech__Fermion *p_fermion;
  double d_result = 0, d_y;

  if( d_ye < 0. || d_ye > 1. )
  {
    fprintf( stderr, "Invalid Ye" );
    exit( EXIT_FAILURE );
  }

  if( strcmp( Libnucnet__Species__getName( p_species ), "n" ) == 0 )
  {

    p_fermion =
      Libstatmech__Fermion__new(
        S_NUCEQ_NEUTRON,
        GSL_CONST_CGSM_UNIFIED_ATOMIC_MASS *
          gsl_pow_2( GSL_CONST_CGSM_SPEED_OF_LIGHT ) *
          GSL_CONST_NUM_MICRO / GSL_CONST_CGSM_ELECTRON_VOLT *
          Libnucnet__Species__getA( p_species ) +
          Libnucnet__Species__getMassExcess( p_species ),
        2,
        0.
      );

    d_y =
      Libstatmech__Fermion__computeQuantity(
        p_fermion,
        S_NUMBER_DENSITY,
        d_T9 * GSL_CONST_NUM_GIGA,
        Libnuceq__getMunkT( p_equil ),
        NULL,
        NULL
      ) / ( d_rho * GSL_CONST_NUM_AVOGADRO );

    d_result =
      log( d_y ) -
      Libnuceq__getMunkT( p_equil ) -
      nse_exp_factor(
        Libnuceq__getNuc( p_equil ),
        p_species,
        d_T9,
        d_rho
      );

    Libstatmech__Fermion__free( p_fermion ) ;

  }
  else if( strcmp( Libnucnet__Species__getName( p_species ), "h1" ) == 0 )
  {

    p_fermion =
      Libstatmech__Fermion__new(
        S_NUCEQ_PROTON,
        GSL_CONST_CGSM_UNIFIED_ATOMIC_MASS *
          gsl_pow_2( GSL_CONST_CGSM_SPEED_OF_LIGHT ) *
          GSL_CONST_NUM_MICRO / GSL_CONST_CGSM_ELECTRON_VOLT *
          Libnucnet__Species__getA( p_species ) +
          Libnucnet__Species__getMassExcess( p_species ),
        2,
        0.
      );

    d_y =
      Libstatmech__Fermion__computeQuantity(
        p_fermion,
        S_NUMBER_DENSITY,
        d_T9 * GSL_CONST_NUM_GIGA,
        Libnuceq__getMupkT( p_equil ),
        NULL,
        NULL
      ) / ( d_rho * GSL_CONST_NUM_AVOGADRO );

    d_result =
      log( d_y ) -
      Libnuceq__getMupkT( p_equil ) -
      nse_exp_factor(
        Libnuceq__getNuc( p_equil ),
        p_species,
        d_T9,
        d_rho
      );

    Libstatmech__Fermion__free( p_fermion ) ;

  }


  return d_result;

}

/*##############################################################################
// nse_exp_factor().
//############################################################################*/

double
nse_exp_factor(
  Libnucnet__Nuc *p_nuc,
  Libnucnet__Species *p_species,
  double d_T9,
  double d_rho
)
{

  return
    log(
      Libnucnet__Species__computeQuantumAbundance(
        p_species,
        d_T9,
        d_rho
      )
    ) -
    Libnucnet__Nuc__computeSpeciesBindingEnergy(
      p_nuc,
      p_species
    ) /
    (
      d_T9 *
      GSL_CONST_NUM_GIGA *
      ( GSL_CONST_NUM_MICRO / GSL_CONST_CGSM_ELECTRON_VOLT ) *
      GSL_CONST_CGSM_BOLTZMANN
    );

}
