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
//       store bosons, compute thermodynamic quantities, and free
//       allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libstatmech.h>
#include "boson_ground_state_functions.h"

/*##############################################################################
// Typedef.
//############################################################################*/

int main( int argc, char **argv )
{

  Libstatmech__Boson * p_boson;
  double d_mukT, d_number_density, d_E1, d_volume;

  if( argc != 7 ) {
    fprintf(
      stderr,
      "\nUsage: %s name rest_mass multiplicity charge temperature number\n\n",
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
      "  temperature = temperature in K\n\n"
    );
    fprintf(
      stderr,
      "  number = number of particles in 1 cc volume\n\n"
    );
    return EXIT_FAILURE;
  } 

  /*============================================================================
  // Set volume in cubic cm.
  //==========================================================================*/

  d_volume = 1.;

  /*============================================================================
  // Calculate number density. 
  //==========================================================================*/

  d_number_density = atof( argv[6] ) / d_volume;

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
  // Update quantities.  Set integral lower limit.
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
      S_PRESSURE,
      (Libstatmech__Boson__Function) boson_pressure_function,
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

  /*============================================================================
  // Compute chemical potential and number density.
  //==========================================================================*/

  d_mukT =
    Libstatmech__Boson__computeChemicalPotential(
      p_boson,
      atof( argv[5] ), 
      d_number_density,
      &d_volume,
      NULL
    );

  fprintf(
    stdout,
    "The mu/kT = %g\n",
    d_mukT
  );

  fprintf(
    stdout,
    "The total number density is %e 1/cm^3\n",
    Libstatmech__Boson__computeQuantity(
      p_boson,
      S_NUMBER_DENSITY,
      atof( argv[5] ),
      d_mukT,
      &d_volume,
      NULL
    )
  );

  /*============================================================================
  // Output populations of ground state level.
  //==========================================================================*/

  fprintf(
    stdout,
    "The number density in ground state is  %e\n",
    boson_number_density_function( 
      p_boson,
      atof( argv[5] ),
      d_mukT,
      &d_volume
    )
  );

  /*============================================================================
  // Output other thermodynamic quantities.
  //==========================================================================*/

  fprintf(
    stdout,
    "The energy density is %e ergs/cm^3\n",
    Libstatmech__Boson__computeQuantity(
      p_boson,
      S_ENERGY_DENSITY,
      atof( argv[5] ),
      d_mukT,
      &d_volume,
      NULL
    )
  );

  fprintf(
    stdout,
    "The internal energy density is %e ergs/cm^3\n",
    Libstatmech__Boson__computeQuantity(
      p_boson,
      S_INTERNAL_ENERGY_DENSITY,
      atof( argv[5] ),
      d_mukT,
      &d_volume,
      NULL
    )
  );

  fprintf(
    stdout,
    "The pressure is %e dynes/cm^2\n",
    Libstatmech__Boson__computeQuantity(
      p_boson,
      S_PRESSURE,
      atof( argv[5] ),
      d_mukT,
      &d_volume,
      NULL
    )
  );

  fprintf(
    stdout,
    "The entropy density is %e ergs/(K*cm^3)\n",
    Libstatmech__Boson__computeQuantity(
      p_boson,
      S_ENTROPY_DENSITY,
      atof( argv[5] ),
      d_mukT,
      &d_volume,
      NULL
    )
  );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  Libstatmech__Boson__free( p_boson );

  return EXIT_SUCCESS;

}
