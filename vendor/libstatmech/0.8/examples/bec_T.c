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
//       store bosons, print out the abundance of the ground state at a 
//       variety of number densities, and free allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libstatmech.h>
#include "boson_ground_state_functions.h"

#define I_STEPS 220
#define D_NMAX 30.         /* maximum number density in powers of 10 */

int main( int argc, char **argv )
{

  Libstatmech__Boson * p_boson;
  double d_n, d_n0, d_mukT, d_E1, d_volume;
  size_t i;

  if( argc != 6 ) {
    fprintf(
      stderr,
      "\nUsage: %s name rest_mass multiplicity charge temperature\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  name = string giving name of the boson\n\n"
    );
    fprintf(
      stderr,
      "  rest_mass = rest mass of boson in MeV\n\n"
    );
    fprintf(
      stderr,
      "  multiplicity = multiplicity of boson\n\n"
    );
    fprintf(
      stderr,
      "  charge = charge of boson\n\n"
    );
    fprintf(
      stderr,
      "  temperature = temperature in K \n\n"
    );
    return EXIT_FAILURE;
  } 

  /*============================================================================
  // Set volume in cubic cm.
  //==========================================================================*/

  d_volume = 1;

  /*============================================================================
  // Create boson.
  //==========================================================================*/

  p_boson =
    Libstatmech__Boson__new(
      argv[1],
      atof( argv[2] ),
      atoi( argv[3] ),
      atof( argv[4] )
   );

  /*============================================================================
  // Compute d_E1. 
  //==========================================================================*/

  d_E1 = 
    2 *
    gsl_pow_2( 
      GSL_CONST_CGSM_PLANCKS_CONSTANT_HBAR *
      M_PI
    ) /
    ( 
      (
        Libstatmech__Boson__getRestMass( p_boson ) *     
        GSL_CONST_CGSM_ELECTRON_VOLT * GSL_CONST_NUM_MEGA / 
        gsl_pow_2( GSL_CONST_CGSM_SPEED_OF_LIGHT )
      ) 
      *
      pow( d_volume, 2./3. )
    );

  /*============================================================================
  // Update quantities.
  //==========================================================================*/

  if(
    !Libstatmech__Boson__updateQuantity(
      p_boson,
      S_NUMBER_DENSITY,
      (Libstatmech__Boson__Function) boson_number_density_function,
      DEFAULT_INTEGRAND
    )
  )
  {
    fprintf( stderr, "Couldn't update function!\n" );
    exit( EXIT_FAILURE );
  }

  if(
    !Libstatmech__Boson__updateQuantity(
      p_boson,
      S_ENTROPY_DENSITY,
      (Libstatmech__Boson__Function) boson_entropy_density_function,
      DEFAULT_INTEGRAND
    )
  )
  {
    fprintf( stderr, "Couldn't update function!\n" );
    exit( EXIT_FAILURE );
  }

  if(
    !Libstatmech__Boson__updateQuantity(
      p_boson,
      S_PRESSURE,
      (Libstatmech__Boson__Function) boson_pressure_function,
      DEFAULT_INTEGRAND
    )
  )
  {
    fprintf( stderr, "Couldn't update function!\n" );
    exit( EXIT_FAILURE );
  }

  /*============================================================================
  // Set integral lower limits.
  //==========================================================================*/

  if(
    !Libstatmech__Boson__updateIntegralLowerLimit(
      p_boson,
      S_NUMBER_DENSITY,
      d_E1 / ( GSL_CONST_CGSM_BOLTZMANN * atof( argv[5] ) )
    )
  )
  {
    fprintf( stderr, "Couldn't update integral lower limit!\n" );
    exit( EXIT_FAILURE );
  }

  if(
    !Libstatmech__Boson__updateIntegralLowerLimit(
      p_boson,
      S_ENTROPY_DENSITY,
      d_E1 / ( GSL_CONST_CGSM_BOLTZMANN * atof( argv[5] ) )
    )
  )
  {
    fprintf( stderr, "Couldn't update integral lower limit!\n" );
    exit( EXIT_FAILURE );
  }

  if(
    !Libstatmech__Boson__updateIntegralLowerLimit(
      p_boson,
      S_PRESSURE,
      d_E1 / ( GSL_CONST_CGSM_BOLTZMANN * atof( argv[5] ) )
    )
  )
  {
    fprintf( stderr, "Couldn't update integral lower limit!\n" );
    exit( EXIT_FAILURE );
  }

  /*============================================================================
  // Output boson info and header.
  //==========================================================================*/

  fprintf(
    stdout,
    "\nFor %s, mc^2 = %g MeV and g = %d and T = %e K\n\n",
    Libstatmech__Boson__getName( p_boson ),
    Libstatmech__Boson__getRestMass( p_boson ),
    Libstatmech__Boson__getMultiplicity( p_boson ),
    atof( argv[5] )
  );


  fprintf(
    stdout,
    "  n (per cc)        N_0/N          mu/kT       E (ergs/cc)    P (dynes/cm^2)   s (ergs / K cm^3)\n"
  );

  /*============================================================================
  // Compute chemical potential and various quantities, and output.
  //==========================================================================*/

  for( i = 0; i <= I_STEPS; i++ )
  {

    d_n = 
      pow( 10., (double) i * D_NMAX / (double) I_STEPS ); 

    d_mukT =
      Libstatmech__Boson__computeChemicalPotential(
        p_boson,
        atof( argv[5] ), 
        d_n, 
        &d_volume,
        NULL
      );

    d_n0 =
      boson_number_density_function( 
        p_boson,
        atof( argv[5] ), 
        d_mukT,
        &d_volume
      );

    fprintf(
      stdout,
      "%12.4e%15.4e%15.4e%15.4e%15.4e%17.4e\n",
      d_n,
      d_n0 / d_n,
      d_mukT,
      Libstatmech__Boson__computeQuantity(
        p_boson,
        S_ENERGY_DENSITY,
        atof( argv[5] ),
        d_mukT,
        &d_volume,
        NULL
      ),
      Libstatmech__Boson__computeQuantity(
        p_boson,
        S_PRESSURE,
        atof( argv[5] ),
        d_mukT,
        &d_volume,
        NULL
      ),
      Libstatmech__Boson__computeQuantity(
        p_boson,
        S_ENTROPY_DENSITY,
        atof( argv[5] ),
        d_mukT,
        &d_volume,
        NULL
      )
    );

  } 

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  Libstatmech__Boson__free( p_boson );

  return EXIT_SUCCESS;

}
