/*//////////////////////////////////////////////////////////////////////////////
// <file type="public">
//   <license>
//     This file was originally written by Tianhong Yu and Bradley S. Meyer.
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
//       Example to demonstrate how to add a Yukawa potential to the integrand, 
//       compute number density, and free
//       allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libstatmech.h>
#include <gsl/gsl_integration.h>

/*##############################################################################
// Global work structure.
//############################################################################*/

typedef struct
{
  double dRadius;
  double dZ;
  double dr_0;
} work;

/*##############################################################################
// Prototypes.
//############################################################################*/

double
n_vs_r_and_x(
  Libstatmech__Fermion *,
  double,
  double,
  double,
  void *
);


double
my_number_density_integrand(
  Libstatmech__Fermion *,
  double,
  double,
  double,
  work *
);

double
potential_function_Yukawa( double, double, double, double );

/*##############################################################################
// main().
//############################################################################*/

int main( int argc, char **argv )
{

  double d_r, d_n_av, d_n_av_check, d_mukT, d_ephi_kT;
  Libstatmech__Fermion * p_fermion, * p_new_fermion;
  work * p_work;

  if( argc != 5 ) {
    fprintf(
      stderr,
      "\nUsage: %s temperature atomic_number atomic_radius r_0\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "temperature = temperature in K\n\n"
    );
    fprintf(
      stderr,
      "atomic_number = atomic number\n\n"
    );
    fprintf(
      stderr,
      "atomic_radius = atomic radius in Bohr radii\n\n"
    );
    fprintf(
      stderr,
      "r_0 = efolding scale of Yukawa potential in Bohr radii\n\n"
    );
    return EXIT_FAILURE;
  } 

  /*============================================================================
  // Create work structure.
  //==========================================================================*/

  p_work = ( work * ) malloc( sizeof( work ) );

  if( !p_work )
  {
    fprintf( stderr, "Couldn't allocate work structure.\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Create fermion. 
  //==========================================================================*/

  p_fermion =
    Libstatmech__Fermion__new(
      "electron", 
      0.511, 
      2, 
      -1.
   );

  /*============================================================================
  // Assign data and print out. 
  //==========================================================================*/

  p_work->dZ= atof( argv[2] );
  p_work->dRadius = atof( argv[3] ) * GSL_CONST_CGSM_BOHR_RADIUS;
  p_work->dr_0 = atof( argv[4] ) * GSL_CONST_CGSM_BOHR_RADIUS;

  /*============================================================================
  // Update integrand and upper limit.
  //==========================================================================*/

  if(
    !Libstatmech__Fermion__updateQuantity(
      p_fermion,
      S_NUMBER_DENSITY,
      NULL,
      (Libstatmech__Fermion__Integrand) my_number_density_integrand
    )
  )
  {
    fprintf( stderr, "Couldn't update function!\n" );
    exit( EXIT_FAILURE );
  }

  Libstatmech__Fermion__updateIntegralLowerLimit(
    p_fermion,
    S_NUMBER_DENSITY,
    1.e-13
  );

  Libstatmech__Fermion__updateIntegralUpperLimit(
    p_fermion,
    S_NUMBER_DENSITY,
    p_work->dRadius
  );

  /*============================================================================
  // Compute chemical potential.
  //==========================================================================*/

  d_n_av = 
    3. * p_work->dZ /
    ( 4. *
      M_PI *
      gsl_pow_3( p_work->dRadius )
    );

  d_mukT =
    Libstatmech__Fermion__computeChemicalPotential(
      p_fermion,
      atof( argv[1] ),
      d_n_av,
      NULL,
      p_work
  );

  d_n_av_check =
    Libstatmech__Fermion__computeQuantity(
      p_fermion,
      S_NUMBER_DENSITY,
      atof( argv[1] ),
      d_mukT,
      NULL,
      p_work
    );

  /*============================================================================
  // Print diagnostics.
  //==========================================================================*/

  fprintf(
    stdout,
    "\nInput:\n\n  T = %s K\n  Z = %s\n  R = %s\n  r_0 = %s\n",
    argv[1],
    argv[2],
    argv[3],
    argv[4]
  );

  fprintf(
    stdout,
    "\nn_av = %g (per cc)\nn_av_check = %g (per cc)\nmu/kT = %g\n\n",
    d_n_av, d_n_av_check,
    d_mukT
  );

  /*============================================================================
  // Print n(r).
  //==========================================================================*/

  p_new_fermion = Libstatmech__Fermion__copy( p_fermion );

  Libstatmech__Fermion__updateQuantity(
    p_new_fermion,
    S_NUMBER_DENSITY,
    NULL,
    (Libstatmech__Fermion__Integrand) n_vs_r_and_x
  );

  fprintf( stdout, "  r/R       ePhi/kT      mu_eff/kT      n (per cc)\n\n" );

  for(
    d_r = 0.02;
    d_r <= 1.01;
    d_r += 0.02
  )
  {

    d_ephi_kT =
      potential_function_Yukawa( 
        atof( argv[2] ),
        atof( argv[1] ),
        d_r * p_work->dRadius,
        p_work->dr_0
      );


    fprintf(
      stdout,
      "%.4f  %.6e  %.6e  %.10e\n",
      d_r,
      d_ephi_kT,
      d_mukT + d_ephi_kT,
      Libstatmech__Fermion__computeQuantity(
        p_new_fermion,
        S_NUMBER_DENSITY,
        atof( argv[1] ),
        d_mukT + d_ephi_kT,
        NULL,
        NULL
      )
    );

  }

  Libstatmech__Fermion__free( p_new_fermion );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  free( p_work );
  Libstatmech__Fermion__free( p_fermion );
  return EXIT_SUCCESS;

}

/*##############################################################################
// n_vs_r_and_x().
//############################################################################*/

double
n_vs_r_and_x(
  Libstatmech__Fermion * p_fermion,
  double d_x,
  double d_T,
  double d_mu_eff_kT,
  void * p_data
)
{

  double d_f, d_gamma, d_mc2;

  if( p_data )
  {
    fprintf( stderr, "No extra data to this routine!\n" );
    exit( EXIT_FAILURE );
  }

  d_mc2 =
    Libstatmech__Fermion__getRestMass( p_fermion ) *
    GSL_CONST_CGSM_ELECTRON_VOLT *
    GSL_CONST_NUM_MEGA;

  d_gamma =
    d_mc2 / ( GSL_CONST_CGSM_BOLTZMANN * d_T );

  if( d_x - d_mu_eff_kT >= 0 )
    d_f =
      exp( d_mu_eff_kT - d_x ) /
      ( exp( d_mu_eff_kT - d_x ) + 1. );
  else
    d_f =
      1. / ( exp( d_x - d_mu_eff_kT ) + 1. );

  /*----------------------------------------------------------------------------
  // Compute integrand.  Check for case that -alpha nearly = gamma.
  //--------------------------------------------------------------------------*/

  d_f *=
    sqrt ( d_x * d_x + 2 * d_x * d_gamma ) *
    ( d_x + d_gamma );

  /*----------------------------------------------------------------------------
  // Apply factor in front.
  //--------------------------------------------------------------------------*/

  d_f *=
    gsl_pow_3( d_mc2 ) *
    Libstatmech__Fermion__getMultiplicity( p_fermion ) /
    (
      2. *
      gsl_pow_2( M_PI ) *
      gsl_pow_3(
        GSL_CONST_CGSM_PLANCKS_CONSTANT_HBAR *
        GSL_CONST_CGSM_SPEED_OF_LIGHT *
        d_gamma
      )
    );

  return d_f;

}

/*##############################################################################
// potential_function_Yukawa() 
//############################################################################*/

double
potential_function_Yukawa(
  double d_Z,
  double d_T,
  double d_r,
  double d_r0
){

  return
    d_Z * 
    GSL_CONST_NUM_FINE_STRUCTURE *
    GSL_CONST_CGSM_PLANCKS_CONSTANT_HBAR *
    GSL_CONST_CGSM_SPEED_OF_LIGHT *
    exp( -d_r /  d_r0 ) 
    /
    (
      d_r *
      GSL_CONST_CGSM_BOLTZMANN *
      d_T
    );

} 

/*##############################################################################
// my_number_density_integrand().
//############################################################################*/

double
my_number_density_integrand(
  Libstatmech__Fermion * p_fermion,
  double d_r,
  double d_T,
  double d_mukT,
  work * p_work
)
{

  Libstatmech__Fermion * p_new_fermion;

  double d_result, d_mueff_kT;

  p_new_fermion = Libstatmech__Fermion__copy( p_fermion );

  d_mueff_kT =
    d_mukT +
    potential_function_Yukawa(
      p_work->dZ,
      d_T,
      d_r,
      p_work->dr_0
    );

  Libstatmech__Fermion__updateQuantity(
    p_new_fermion,
    S_NUMBER_DENSITY,
    NULL,
    (Libstatmech__Fermion__Integrand) n_vs_r_and_x
  );

  /*============================================================================
  // Uncomment the lines below if calculation can't converge.
  //==========================================================================*/
  
/*
  Libstatmech__Fermion__updateQuantityIntegralAccuracy(
    p_new_fermion,
    S_NUMBER_DENSITY,
    1.e-8,
    1.e-2
  );
*/

  d_result =
    3. *
    gsl_pow_2( d_r ) *
    Libstatmech__Fermion__computeQuantity(
      p_new_fermion,
      S_NUMBER_DENSITY,
      d_T,
      d_mueff_kT,
      NULL,
      NULL
    )
    / gsl_pow_3( p_work->dRadius );

  Libstatmech__Fermion__free( p_new_fermion );

  return d_result;

}

