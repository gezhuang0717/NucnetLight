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
//       Example to demonstrate how to use libstatmech routines to create and
//       store fermions, compute thermodynamic quantities, and free
//       allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libstatmech.h>

/*##############################################################################
// Prototypes. 
//############################################################################*/

double
my_pressure_integrand(
  Libstatmech__Fermion *, double, double, double, void *
);

double
my_pressure_function(
  Libstatmech__Fermion *, double, double, void *
);

int main( int argc, char **argv )
{

  Libstatmech__Fermion * p_fermion;
  double d_mukT;

  if( argc != 7 ) {
    fprintf(
      stderr,
      "\nUsage: %s name rest_mass multiplicity charge temperature number_density\n\n",
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
      "  number_density = number density in per cc\n\n"
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
  // Compute chemical potential and print it out. 
  //==========================================================================*/

  d_mukT =
    Libstatmech__Fermion__computeChemicalPotential(
      p_fermion,
      atof( argv[5] ),
      atof( argv[6] ),
      NULL,
      NULL
    );

  fprintf(
    stdout,
    "The mu/kT = %g\n",
    d_mukT
  );

  /*============================================================================
  // Print other thermal quantities. 
  //==========================================================================*/

  fprintf(
    stdout,
    "The number density is %22.16e cm^-3\n",
    Libstatmech__Fermion__computeQuantity(
      p_fermion,
      S_NUMBER_DENSITY,
      atof( argv[5] ),
      d_mukT,
      NULL, 
      NULL
    )
  );

  fprintf(
    stdout,
    "The pressure is %22.16e dynes/cm^2\n",
    Libstatmech__Fermion__computeQuantity(
      p_fermion,
      S_PRESSURE,
      atof( argv[5] ),
      d_mukT,
      NULL, 
      NULL
    )
  );

  fprintf(
    stdout,
    "The energy density is %e ergs/cm^3\n",
    Libstatmech__Fermion__computeQuantity(
      p_fermion,
      S_ENERGY_DENSITY,
      atof( argv[5] ),
      d_mukT,
      NULL,
      NULL
    )
  );

  fprintf(
    stdout,
    "The internal energy density is %e ergs/cm^3\n",
    Libstatmech__Fermion__computeQuantity(
      p_fermion,
      S_INTERNAL_ENERGY_DENSITY,
      atof( argv[5] ),
      d_mukT,
      NULL,
      NULL
    )
  );

  fprintf(
    stdout,
    "The entropy density is %e ergs/(cm^3*K)\n",
    Libstatmech__Fermion__computeQuantity(
      p_fermion,
      S_ENTROPY_DENSITY,
      atof( argv[5] ),
      d_mukT,
      NULL,
      NULL
    )
  );

  /*============================================================================
  // Update pressure integrand and print the new pressure.
  //==========================================================================*/

  fprintf( stdout, "\nCompute the pressure by updating the quantity:\n" );

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

  fprintf(
    stdout,
    "With our integrand the pressure is %e\n",
    Libstatmech__Fermion__computeQuantity(
      p_fermion,
      S_PRESSURE,
      atof( argv[5] ),
      d_mukT,
      NULL,
      NULL
    )
  );

  /*============================================================================
  // Compute the pressure with the user-defined function.
  //==========================================================================*/

  if(
    !Libstatmech__Fermion__updateQuantity(
      p_fermion,
      S_PRESSURE,
      (Libstatmech__Fermion__Function) my_pressure_function,
      NULL
    )
  )
  {
    fprintf( stderr, "Couldn't update function!\n" );
    exit( EXIT_FAILURE );
  }

  fprintf(
    stdout,
    "With our function the pressure is %e\n",
    Libstatmech__Fermion__computeQuantity(
      p_fermion,
      S_PRESSURE,
      atof( argv[5] ),
      d_mukT,
      NULL,
      NULL
    )
  );

  /*============================================================================
  // Repeat the last two calculations by instead using user-defined quantities.
  //==========================================================================*/

  fprintf( stdout, "\nCompute the pressure with user-defined quantities:\n" );

  if( 
      !Libstatmech__Fermion__updateQuantity(
         p_fermion,
         "my pressure integrand",
         NULL,
         (Libstatmech__Fermion__Integrand) my_pressure_integrand
      )
  )
  {
     fprintf( stderr, "Couldn't create quantity.\n" );
     return EXIT_FAILURE;
  }

  if( 
      !Libstatmech__Fermion__updateQuantity(
         p_fermion,
         "my pressure function",
         (Libstatmech__Fermion__Function) my_pressure_function,
         NULL
      )
  )
  {
     fprintf( stderr, "Couldn't create quantity.\n" );
     return EXIT_FAILURE;
  }

  fprintf(
    stdout,
    "With our integrand the pressure is %e\n",
    Libstatmech__Fermion__computeQuantity(
      p_fermion,
      "my pressure integrand",
      atof( argv[5] ),
      d_mukT,
      NULL,
      NULL
    )
  );

  fprintf(
    stdout,
    "With our function the pressure is %e\n",
    Libstatmech__Fermion__computeQuantity(
      p_fermion,
      "my pressure function",
      atof( argv[5] ),
      d_mukT,
      NULL,
      NULL
    )
  );

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

/*##############################################################################
// my_pressure_function()
//############################################################################*/

double
my_pressure_function(
  Libstatmech__Fermion *p_fermion,
  double d_T,
  double d_alpha,
  void *p_user_data
){

  double d_result;

  if( p_user_data ) {
    fprintf(
      stderr,
      "No user data should be passed to this routine!\n"
    );
    exit( EXIT_FAILURE );
  }

  d_result = 
    Libstatmech__Fermion__computeQuantity(
      p_fermion,
      S_NUMBER_DENSITY,
      d_T,
      d_alpha,
      NULL,
      NULL
    ) * GSL_CONST_CGSM_BOLTZMANN * d_T ;

  if ( !gsl_finite( d_result ) )
  {
    fprintf( stderr, "Infinite function.\n" ); 
    exit( EXIT_FAILURE );
  } 

  return d_result;

}

