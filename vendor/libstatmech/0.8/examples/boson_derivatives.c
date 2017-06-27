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
//       store bosons, compute temperature derivatives of thermodynamic
//       quantities, and free allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libstatmech.h>

int main( int argc, char **argv )
{

  Libstatmech__Boson * p_boson;
  double d_mukT;

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
      "  number = number of particles\n\n"
    );
    return EXIT_FAILURE;
  } 

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
  // Compute chemical potential and output.
  //==========================================================================*/

  d_mukT =
    Libstatmech__Boson__computeChemicalPotential(
      p_boson, atof( argv[5] ), atof( argv[6] ), NULL, NULL
    );

  fprintf(
    stdout,
    "The mu/kT = %g\n",
    d_mukT
  );

  /*============================================================================
  // Compute temperature derivatives and heat capacity. 
  //==========================================================================*/

  fprintf(
    stdout,
    "The d d_mu_kT / d d_T is %e K^-1 \n",
    Libstatmech__Boson__computeTemperatureDerivative(
      p_boson,
      S_CHEMICAL_POTENTIAL,
      atof( argv[5] ),
      atof( argv[6] ),
      NULL,
      NULL 
    )
  );

  fprintf(
    stdout,
    "The temperature derivative of pressure is %e dynes/(cm^2*K)\n",
    Libstatmech__Boson__computeTemperatureDerivative(
      p_boson,
      S_PRESSURE,
      atof( argv[5] ),
      atof( argv[6] ),
      NULL,
      NULL 
    )
  );

  fprintf(
    stdout,
    "The temperature derivative of energy density is %e ergs/(cm^3*K)\n",
    Libstatmech__Boson__computeTemperatureDerivative(
      p_boson,
      S_ENERGY_DENSITY,
      atof( argv[5] ),
      atof( argv[6] ),
      NULL,
      NULL 
    )
  );

  fprintf(
    stdout,
    "The temperature derivative of internal energy density is %e ergs/(cm^3*K)\n",
    Libstatmech__Boson__computeTemperatureDerivative(
      p_boson,
      S_INTERNAL_ENERGY_DENSITY,
      atof( argv[5] ),
      atof( argv[6] ),
      NULL,
      NULL 
    )
  );

  fprintf(
    stdout,
    "The temperature derivative of the entropy density is %e ergs/(cm^3*K^2)\n",
    Libstatmech__Boson__computeTemperatureDerivative(
      p_boson,
      S_ENTROPY_DENSITY,
      atof( argv[5] ),
      atof( argv[6] ),
      NULL,
      NULL
    )
  );

  fprintf(
    stdout,
    "The heat capacity density is %e ergs/(cm^3*K)\n",
    atof( argv[5] ) *
    Libstatmech__Boson__computeTemperatureDerivative(
      p_boson,
      S_ENTROPY_DENSITY,
      atof( argv[5] ),
      atof( argv[6] ),
      NULL,
      NULL
    )
  );

  /*============================================================================
  // Clean up and exit. 
  //==========================================================================*/

  Libstatmech__Boson__free( p_boson );
  return EXIT_SUCCESS;

}
