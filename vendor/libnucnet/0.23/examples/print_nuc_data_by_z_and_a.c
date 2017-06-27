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
//       Example to demonstrate how to use libnucnet routines to parse
//       in a Libnucnet__Nuc structure of nuclear species from an input
//       xml file, print out the index number of a species selected by
//       its atomic number, mass number, and state,
//       and clear the structure and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libnucnet__Nuc.h>

int main( int argc, char * argv[] ) {

  Libnucnet__Nuc *p_my_nuclei;
  Libnucnet__Species *p_species;

  /*============================================================================
  // Check input.
  //==========================================================================*/


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
       stderr, "  State = State of species\n\n"
     ); 
     return EXIT_FAILURE;
   }

  /*============================================================================
  // Get network from input xml file.
  //==========================================================================*/

  p_my_nuclei = Libnucnet__Nuc__new_from_xml( argv[1], NULL );

  /*============================================================================
  // Get species index and name, if present.
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
  // Print out data if species present or indicate species not in network.
  //==========================================================================*/

  if( p_species ) {
    fprintf(
      stdout,
      "Index of %s = %lu\n",
      Libnucnet__Species__getName( p_species ),
      (unsigned long) Libnucnet__Species__getIndex( p_species )
    );
    Libnucnet__Nuc__free( p_my_nuclei );
    return EXIT_SUCCESS;
  } else {
    fprintf( stderr, "Species not present in network!\n" );
    Libnucnet__Nuc__free( p_my_nuclei );
    return EXIT_FAILURE;
  }

}
