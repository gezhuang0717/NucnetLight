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
//       Example to demonstrate how to compute the pressure for a gas of
//       fermions with the default and with a user-supplied pressure integrand,
//       and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libstatmech.h>

/*##############################################################################
// Constant definition. 
//############################################################################*/

#define   D_FACTOR   30.

/*##############################################################################
// Prototype. 
//############################################################################*/

double
my_pressure_integrand(
  Libstatmech__Fermion *, double, double, double, void *
);

int main( int argc, char **argv )
{

  Libstatmech__Fermion * p_fermion;
  double d_x, d_default[1000], d_my_pressure[1000];
  int i;

  if( argc != 7 ) {
    fprintf(
      stderr,
      "\nUsage: %s name rest_mass multiplicity charge temperature mukT\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  name = string giving name of the fermion\n\n"
    );
    fprintf(
      stderr,
      "  rest_mass = rest mass of fermion in MeV\n\n"
    );
    fprintf(
      stderr,
      "  multiplicity = multiplicity of fermion\n\n"
    );
    fprintf(
      stderr,
      "  charge = charge of fermion\n\n"
    );
    fprintf(
      stderr,
      "  temperature = temperature in K\n\n"
    );
    fprintf(
      stderr,
      "  mukT = mu/kT for the fermion\n\n"
    );
    return EXIT_FAILURE;
  } 

  /*============================================================================
  // Create fermion. 
  //==========================================================================*/

  p_fermion =
    Libstatmech__Fermion__new(
      argv[1],
      atof( argv[2] ),
      atoi( argv[3] ),
      atof( argv[4] )
   );

  /*============================================================================
  // First use user-defined pressure integrand.
  //==========================================================================*/

  if(
    !Libstatmech__Fermion__updateQuantity(
      p_fermion,
      S_PRESSURE,
      NULL,
      (Libstatmech__Fermion__Integrand) my_pressure_integrand
    )
  )
  {
    fprintf( stderr, "Couldn't update quantity!\n" );
    exit( EXIT_FAILURE );
  }

  d_x = 0.;

  for( i = 0; i < 100; i++ )
  {
    d_my_pressure[i] =
      Libstatmech__Fermion__computeIntegrandValue(
        p_fermion,
        S_PRESSURE,
        d_x,
        atof( argv[5] ),
        atof( argv[6] ),
        NULL
      );
    d_x += fabs( atof( argv[6] ) ) / D_FACTOR;
  }
        
  /*============================================================================
  // Now reset to default pressure.
  //==========================================================================*/

  if(
     !(
        Libstatmech__Fermion__updateQuantity(
          p_fermion,
          S_PRESSURE,
          NULL,
          DEFAULT_INTEGRAND
        )
     )
  )
  {
     fprintf( stderr, "Couldn't reset to default.\n" );
     return EXIT_FAILURE;
  }

  /*============================================================================
  // Now compute default pressure integrand.
  //==========================================================================*/

  d_x = 0.;

  for( i = 0; i < 100; i++ )
  {
    d_default[i] =
      Libstatmech__Fermion__computeIntegrandValue(
        p_fermion,
        S_PRESSURE,
        d_x,
        atof( argv[5] ),
        atof( argv[6] ),
        NULL
      );
    d_x += fabs( atof( argv[6] ) ) / D_FACTOR;
  }
        
  /*============================================================================
  // Now print results.
  //==========================================================================*/

  fprintf(
    stdout,
    "      x         Default P(x)      My P(x)\n"
  );

  d_x = 0.;

  for( i = 0; i < 100; i++ )
  {
    fprintf(
      stdout,
      "%12.4e%15.4e%15.4e\n",
      d_x,
      d_default[i],
      d_my_pressure[i]
    );
    d_x += fabs( atof( argv[6] ) ) / D_FACTOR;
  }

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  Libstatmech__Fermion__free( p_fermion );
  return EXIT_SUCCESS;

}

/*##############################################################################
// my_pressure_integrand()
//############################################################################*/

double
my_pressure_integrand(
  Libstatmech__Fermion *p_fermion,
  double d_x,
  double d_T,
  double d_alpha,
  void *p_user_data
){

  double d_f, d_gamma, d_mc2;

  if( p_user_data ) {
    fprintf(
      stderr,
      "No user data should be passed to this routine!\n"
    );
    exit( EXIT_FAILURE );
  }

  d_mc2 =
    Libstatmech__Fermion__getRestMass( p_fermion ) *
    GSL_CONST_CGSM_ELECTRON_VOLT *
    GSL_CONST_NUM_MEGA;

  d_gamma =
    d_mc2 / ( GSL_CONST_CGSM_BOLTZMANN * d_T );

  d_f =
     d_gamma * sqrt ( 2 * d_x * d_gamma ) * exp( d_alpha - d_x );

  d_f *=
    gsl_pow_4( d_mc2 ) /
    (
      gsl_pow_2( M_PI ) *
      gsl_pow_3(
        GSL_CONST_CGSM_PLANCKS_CONSTANT_HBAR *
        GSL_CONST_CGSM_SPEED_OF_LIGHT
      ) *
      gsl_pow_4( d_gamma )
    );

  if ( !gsl_finite( d_f ) )
  {
    fprintf( stderr, "Infinite integrand.\n" ); 
    exit( EXIT_FAILURE );
  } 

  return d_f;

}
