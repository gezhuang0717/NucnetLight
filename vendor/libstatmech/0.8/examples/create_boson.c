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
//       Example to demonstrate how to create and store a boson, print out
//       data about the boson, and free data and boson.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libstatmech.h>

int main( int argc, char **argv )
{

  Libstatmech__Boson * p_boson;

  if( argc != 5 ) {
    fprintf(
      stderr,
      "\nUsage: %s name rest_mass multiplicity charge\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  name = name of the boson\n\n"
    );
    fprintf(
      stderr,
      "  rest_mass = rest mass of boson in MeV\n\n"
    );
    fprintf(
      stderr,
      "  multiplicity = multiplicity of boson (2*spin + 1)\n\n"
    );
    fprintf(
      stderr,
      "  charge = charge of the boson (in units of fundamental charge)\n\n"
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
  // Print properties of boson. 
  //==========================================================================*/

  printf( 
    "The name of the boson is %s\n",
    Libstatmech__Boson__getName( p_boson )
  );

  printf(
    "The rest mass of the boson is %e\n", 
    Libstatmech__Boson__getRestMass( p_boson )
  ); 

  printf(
    "The multiplicity of the boson is %d\n", 
    Libstatmech__Boson__getMultiplicity( p_boson )
  ); 

  printf(
    "The charge of the boson is %e\n", 
    Libstatmech__Boson__getCharge( p_boson )
  ); 

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  Libstatmech__Boson__free( p_boson );
  return EXIT_SUCCESS;

}

