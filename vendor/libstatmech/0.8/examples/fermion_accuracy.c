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
//       store fermions, compute the chemical potential for a number of
//       relative accuracies, and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libstatmech.h>

int main( int argc, char **argv )
{

  Libstatmech__Fermion * p_fermion;
  double d_eps_rel;
  size_t i;

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
  // Compute chemical potential for a variety of error limits and print it out. 
  //==========================================================================*/

  for( i = 0; i <= 13; i++ )
  {

    d_eps_rel = pow( 10., -(double) i );

    Libstatmech__Fermion__updateQuantityIntegralAccuracy(
      p_fermion,
      S_NUMBER_DENSITY,
      0.,
      d_eps_rel
    );

    fprintf(
      stdout,
      "eps_rel = %.2e  mu/kT = %.15e\n",
      d_eps_rel,
      Libstatmech__Fermion__computeChemicalPotential(
        p_fermion,
        atof( argv[5] ),
        atof( argv[6] ),
        NULL,
        NULL
      )
    );

  }

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  Libstatmech__Fermion__free( p_fermion );
  return EXIT_SUCCESS;

}
