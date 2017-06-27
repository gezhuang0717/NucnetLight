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
//       a vector from XML input, print out the data, and free the
//       vector.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include "WnMatrix.h"

int main( int argc, char *argv[] ) {

  gsl_vector *p_my_vector;

  /*============================================================================
  // Check arguments.
  //==========================================================================*/

  if ( argc != 2 && argc != 3 ) {
    fprintf(
      stderr, "\nUsage: %s filename xpath\n\n", argv[0]
    );
    fprintf(
      stderr, "  filename = input xml data file\n\n"
    );
    fprintf(
      stderr, "  xpath = xpath expression (optional)\n\n"
    );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Validate input XML.
  //==========================================================================*/

  if( !WnMatrix__is_valid_vector_input_xml( argv[1] ) ) {
    fprintf( stderr, "Not valid input XML!\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Read in matrix from XML.
  //==========================================================================*/

  if( argc == 2 ) {
    p_my_vector = WnMatrix__new_gsl_vector_from_xml( argv[1], NULL );
  } else {
    p_my_vector = WnMatrix__new_gsl_vector_from_xml( argv[1], argv[2] );
  }

  /*============================================================================
  // Print out.
  //==========================================================================*/

  fprintf( stdout, "\nValue\n" );
  fprintf( stdout, "-----\n" );

  gsl_vector_fprintf( stdout, p_my_vector, "%4.1f" );

  fprintf( stdout, "\n" );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  gsl_vector_free( p_my_vector );

  return EXIT_SUCCESS;

}
