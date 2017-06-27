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
//       a matrix from XML input, convert to coordinate matrix form,
//       print out the data, and free the matrices.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include "WnMatrix.h"

int main( int argc, char *argv[] ) {

  unsigned int i;
  WnMatrix * p_my_matrix;
  WnMatrix__Coo *p_coo;

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

  if( !WnMatrix__is_valid_input_xml( argv[1] ) ) {
    fprintf( stderr, "Not valid input XML!\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Read in matrix from XML.
  //==========================================================================*/

  if( argc == 2 ) {
    p_my_matrix = WnMatrix__new_from_xml( argv[1], NULL );
  } else {
    p_my_matrix = WnMatrix__new_from_xml( argv[1], argv[2] );
  }

  /*============================================================================
  // Get coordinate matrix.
  //==========================================================================*/

  p_coo = WnMatrix__getCoo( p_my_matrix );

  /*============================================================================
  // Print matrix in coordinate format.
  //==========================================================================*/

  fprintf( stdout, "Row   Column   Value\n" );
  fprintf( stdout, "---   ------   -----\n" );

  for( i = 0; i < WnMatrix__getNumberOfElements( p_my_matrix ); i++ ) {

     fprintf(
       stdout,
       "%3lu %8lu %7.1f\n", 
       (unsigned long) ( WnMatrix__Coo__getRowVector( p_coo ) )[i],
       (unsigned long) ( WnMatrix__Coo__getColumnVector( p_coo ) )[i],
       ( WnMatrix__Coo__getValueVector( p_coo ) )[i]
     );

  }

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  WnMatrix__Coo__free( p_coo );

  WnMatrix__free( p_my_matrix );

  return EXIT_SUCCESS;

}
