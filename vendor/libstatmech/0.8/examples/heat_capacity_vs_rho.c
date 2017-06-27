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
//       store fermions, compute the heat capacity at a variety of densities,
//       and free allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libstatmech.h>

int main( int argc, char **argv )
{

  Libstatmech__Fermion * p_fermion;
  double d_mukT, d_cv, d_density, d_number_density;
  int i;

  if( argc != 7 ) {
    fprintf(
      stderr,
      "\nUsage: %s name rest_mass multiplicity charge temperature ye\n\n",
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
      "  ye = electron-to-baryon ratio\n\n"
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
  // Compute heat capacity at a variety of densities. 
  //==========================================================================*/

  fprintf(
    stdout,
    "     rho (g/cc)          mu/kT         C_V (ergs/K*gram)\n"
  );

  for( i = 2; i < 17; i++ )
  {

    d_density = pow( 10., (double) i );

    d_number_density = GSL_CONST_NUM_AVOGADRO * atof( argv[6] ) * d_density;

    d_mukT =
      Libstatmech__Fermion__computeChemicalPotential(
        p_fermion,
        atof( argv[5] ),
        d_number_density,
        NULL,
        NULL
      );

    d_cv =
      atof( argv[5] ) *
      Libstatmech__Fermion__computeTemperatureDerivative(
        p_fermion,
        S_ENTROPY_DENSITY,
        atof( argv[5] ),
        d_number_density,
        NULL,
        NULL
     ) / d_density;

    fprintf( stdout, "%15.4e  %15.4e  %15.4e\n", d_density, d_mukT, d_cv );

  }

  /*============================================================================
  // Clean up and exit. 
  //==========================================================================*/

  Libstatmech__Fermion__free( p_fermion );
  return EXIT_SUCCESS;

}
