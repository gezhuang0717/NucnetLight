/*//////////////////////////////////////////////////////////////////////////////
// <file type="public">
//   <license>
//     This file was originally written by Bradley S. Meyer.
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
//       Example to demonstrate how to use libnucnet routines to create a
//       Libnucnet__Nuc structure of nuclear species from an input xml
//       file, print out the nuclear partition function for a species
//       selected by its Z, A, and state, and clear the structure and free the
//       allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libnucnet__Nuc.h>

int main( int argc, char * argv[] ) {

  Libnucnet__Nuc *p_my_nuclei;
  Libnucnet__Species *p_species;
  size_t i;
  double d_t9;

  /* Check input */

  if ( (argc < 4) | (argc > 5) ) {
       fprintf(
        stderr, "\nUsage: %s filename Z A State\n\n", argv[0]
     );
     fprintf(
       stderr, "  filename = input nuclear data xml filename or url\n\n"
     ); 
     fprintf(
       stderr, "  Z = Z of species\n\n"
     ); 
     fprintf(
       stderr, "  A = A of species\n\n"
     ); 
     fprintf(
       stderr, "  State = State of species (optional)\n\n"
     ); 
     return EXIT_FAILURE;
   }

  /*============================================================================
  // Get network from input xml file.
  //==========================================================================*/

  p_my_nuclei = Libnucnet__Nuc__new_from_xml( argv[1], NULL );

  /*============================================================================
  // Get species index and name.
  //==========================================================================*/

  if( !argv[4] ) {
    p_species =
      Libnucnet__Nuc__getSpeciesByZA(
        p_my_nuclei,
        (unsigned int) atoi( argv[2] ),
        (unsigned int) atoi( argv[3] ),
        NULL
      );
  } else {
    p_species =
      Libnucnet__Nuc__getSpeciesByZA(
        p_my_nuclei,
        (unsigned int) atoi( argv[2] ),
        (unsigned int) atoi( argv[3] ),
        argv[4]
    );
  }

  /*============================================================================
  // Exit if species not present.
  //==========================================================================*/

  if( !p_species ) {
    fprintf( stderr, "Species not present in network!\n" );
    Libnucnet__Nuc__free( p_my_nuclei );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Print out partition function if species present in network.
  //==========================================================================*/

  printf( "\nFor %s\n\n", Libnucnet__Species__getName( p_species ) );

  printf( "      t9        Partition Function\n" );
  printf( "    ------      ------------------\n\n" );

  for( i = 0; i < 31; i++ ) {
    d_t9 = pow( 10., -1 + ( (float) i ) /10. );
    printf("  %8.4f    %15.4e\n", d_t9,
      Libnucnet__Species__computePartitionFunction(
        p_species, d_t9 )
    );
  }

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  Libnucnet__Nuc__free( p_my_nuclei );

  return EXIT_SUCCESS;

}
