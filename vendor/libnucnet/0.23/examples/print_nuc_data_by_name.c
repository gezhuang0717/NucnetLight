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
//       xml file, print out the data about a particular species
//       selected by its name,
//       and clear the structure and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libnucnet__Nuc.h>

int main( int argc, char * argv[] ) {

  size_t i;
  Libnucnet__Nuc *p_my_nuclei;
  Libnucnet__Species *p_species;
  gsl_vector *p_t9 = NULL, *p_log10_partf;

  if ( argc != 3 ) {
       fprintf(
        stderr, "\nUsage: %s filename nucname\n\n", argv[0]
     );
     fprintf(
       stderr, "  filename = input nuclear data xml filename or url\n\n"
     ); 
     fprintf(
       stderr, "  nucname = name of species\n\n"
     ); 
     return EXIT_FAILURE;
   }

  /*============================================================================
  // Get network from input xml file.  For illustration, do this the hard
  // way--create the structure, then update the data.  The alternative,
  // of course, is to call Libnucnet__Nuc__new_from_xml().
  //==========================================================================*/

  p_my_nuclei = Libnucnet__Nuc__new();

  Libnucnet__Nuc__updateFromXml( p_my_nuclei, argv[1], NULL );

  /*============================================================================
  // Check that species present in network.  If not, exit.
  //==========================================================================*/

  if(
    !(
      p_species = Libnucnet__Nuc__getSpeciesByName( p_my_nuclei, argv[2] )
     )
  )
  {
    fprintf( stderr, "%s not in network!\n", argv[2] );
    Libnucnet__Nuc__free( p_my_nuclei );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Print out data if species present in network.
  //==========================================================================*/

  fprintf( stdout, "For %s\n", argv[2] ); 
  fprintf(
    stdout,
    "  Z = %u\n",
    Libnucnet__Species__getZ( p_species )
  );
  fprintf(
    stdout,
    "  A = %u\n",
    Libnucnet__Species__getA( p_species )
  );
  fprintf(
    stdout,
    "  Mass Excess (MeV) = %f\n",
    Libnucnet__Species__getMassExcess( p_species )
  );
  fprintf(
    stdout,
    "  Spin = %f\n",
    Libnucnet__Species__getSpin( p_species )
  );
  if( strlen( Libnucnet__Species__getSource( p_species ) ) != 0 ) {
    printf( "  Data Source: %s\n", Libnucnet__Species__getSource( p_species ) );
  }

  /*============================================================================
  // Partition function data.
  //==========================================================================*/

  p_t9 = Libnucnet__Species__getPartitionFunctionT9( p_species );
  p_log10_partf = Libnucnet__Species__getPartitionFunctionLog10( p_species );

  if( p_t9 )
  {

    fprintf( stdout, "\nPartition function data:\n\n" );

    fprintf( stdout, "Index      T9   Log10 Partition Function / (2J+1)\n" );

    fprintf( stdout, "-----     ----  ---------------------------------\n" );

    for( i = 0; i < WnMatrix__get_gsl_vector_size( p_t9 ); i++ )
      fprintf(
        stdout,
        "%5lu  %7.2f  %12.4e\n",
        (unsigned long) i,
        gsl_vector_get( p_t9, i ),
        gsl_vector_get( p_log10_partf, i )
      );

  }

  /*============================================================================
  // Clean up and return.
  //==========================================================================*/

  Libnucnet__Nuc__free( p_my_nuclei );

  return EXIT_SUCCESS;

}
