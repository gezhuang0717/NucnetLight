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
//       Example to demonstrate how to use libnucnet routines to create
//       a new Libnucnet structure from an input xml file, print out 
//       abundances for the given zone by nucleon number, and clear the
//       structure and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libnucnet.h>

int main( int argc, char * argv[] ) {

  Libnucnet *p_my_nucnet;
  Libnucnet__Zone *p_zone;

  gsl_vector * p_vector;
  size_t i;

  if ( argc != 7 ) {
      fprintf(
        stderr,
        "\nUsage: %s nuc_file mass_file label1 label2 label3 nucleon_type\n\n", argv[0]
      );
      fprintf(
        stderr, "  nuc_file = input nuclear network data xml filename\n\n"
      );
      fprintf(
        stderr, "  mass_file = input mass fraction xml filename\n\n"
      );
      fprintf(
        stderr, "  label1 = first zone label\n\n"
      );
      fprintf(
        stderr, "  label2 = second zone label\n\n"
      );
      fprintf(
        stderr, "  label3 = third zone label\n\n"
      );
      fprintf(
        stderr, "  nucleon_type = nucleon type for print out (z, n, or a)\n\n"
      );

      return EXIT_FAILURE;
  }

  /*============================================================================
  // Read input file.
  //==========================================================================*/

  p_my_nucnet = Libnucnet__new();

  Libnucnet__Nuc__updateFromXml(
    Libnucnet__Net__getNuc(
      Libnucnet__getNet( p_my_nucnet )
    ),
    argv[1],
    NULL
  );

  Libnucnet__assignZoneDataFromXml( p_my_nucnet, argv[2], NULL );

  /*============================================================================
  // Retrieve zone and print abundances.
  //==========================================================================*/

  printf( "\n" );

  p_zone = Libnucnet__getZoneByLabels( 
    p_my_nucnet, argv[3], argv[4], argv[5] 
  );

  if( !p_zone ) {

    fprintf( stderr, "Zone not found!\n" );
    return EXIT_FAILURE;

  }

  p_vector = Libnucnet__Zone__getSummedAbundances( p_zone, argv[6] );

  fprintf(
    stdout,
    "For zone with labels %s %s %s:\n\n",
    Libnucnet__Zone__getLabel( p_zone, 1 ),
    Libnucnet__Zone__getLabel( p_zone, 2 ),
    Libnucnet__Zone__getLabel( p_zone, 3 )
  );

  fprintf( stdout, "%s     Y(%s)\n\n", argv[6], argv[6] );

  for( i = 0; i < WnMatrix__get_gsl_vector_size( p_vector ); i++ )
    fprintf(
      stdout,
      "%lu  %.4e\n",
      (unsigned long) i,
      gsl_vector_get( p_vector, i )
    );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  gsl_vector_free( p_vector );
  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
