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
//     Please see the src/README.txt in this distribution for more copyright
//     and license information.
//   </license>
//   <description>
//     <abstract>
//       Example to demonstrate how to use wn_matrix routines to create
//       and store elements in a vector from an input ascii file, 
//       output the data to an XML file, and free the vector.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include "WnMatrix.h"

int main( int argc, char *argv[] ) {

  gsl_vector *p_my_vector;
  FILE *p_file;
  double d_value;
  size_t i_length = 0;

  /*============================================================================
  // Check arguments.
  //==========================================================================*/

  if ( argc != 3 && argc != 4 ) {
    fprintf(
      stderr, "\nUsage: %s in_file out_file \n\n", argv[0]
    );
    fprintf(
      stderr, "  in_file = input vector txt data file\n\n"
    );
    fprintf(
      stderr, "  out_file = output vector xml data file\n\n"
    );
    fprintf(
      stderr,
      "  format = format code for output of vector element value (optional--%%g if not supplied)\n\n"
    );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Open vector text file.
  //==========================================================================*/

  p_file = fopen( argv[1], "r" );

  if( !p_file ) {
    fprintf(
      stderr,
      "Couldn't open file!\n"
    );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Count lines in vector file.
  //==========================================================================*/

  while( !feof( p_file ) ) {
    fscanf( p_file, "%lf\n", &d_value );
    i_length++;
  }

  /*============================================================================
  // Rewind.
  //==========================================================================*/

  rewind( p_file );

  /*============================================================================
  // Create vector.
  //==========================================================================*/

  p_my_vector = gsl_vector_alloc( i_length );

  /*============================================================================
  // Read in vector data.
  //==========================================================================*/

  gsl_vector_fscanf( p_file, p_my_vector );

  /*============================================================================
  // Close text file.
  //==========================================================================*/

  fclose( p_file );

  /*============================================================================
  // Output vector XML file.
  //==========================================================================*/

  if( argc == 3 )
    WnMatrix__write_gsl_vector_to_xml_file( p_my_vector, argv[2], NULL );
  else
    WnMatrix__write_gsl_vector_to_xml_file( p_my_vector, argv[2], argv[3] );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  gsl_vector_free( p_my_vector );

  return EXIT_SUCCESS;

}
