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
//       Example to demonstrate how to create and store a fermion, print out
//       data about the fermion, and free data and fermion.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libstatmech.h>

int main( int argc, char **argv )
{

  Libstatmech__Fermion * p_fermion;

  if( argc != 5 ) {
    fprintf(
      stderr,
      "\nUsage: %s name rest_mass multiplicity charge\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  name = name of the fermion\n\n"
    );
    fprintf(
      stderr,
      "  rest_mass = rest mass of fermion in MeV\n\n"
    );
    fprintf(
      stderr,
      "  multiplicity = multiplicity of fermion (2*spin + 1)\n\n"
    );
    fprintf(
      stderr,
      "  charge = charge of the fermion (in units of fundamental charge)\n\n"
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
  // Print properties of the fermion. 
  //==========================================================================*/

  printf( 
    "The name of the fermion is %s\n",
    Libstatmech__Fermion__getName( p_fermion )
  );

  printf(
    "The rest mass of the fermion is %e\n", 
    Libstatmech__Fermion__getRestMass( p_fermion )
  ); 

  printf(
    "The multiplicity of the fermion is %d\n", 
    Libstatmech__Fermion__getMultiplicity( p_fermion )
  ); 

  printf(
    "The charge of the fermion is %e\n", 
    Libstatmech__Fermion__getCharge( p_fermion )
  ); 

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  Libstatmech__Fermion__free( p_fermion );
  return EXIT_SUCCESS;

}

