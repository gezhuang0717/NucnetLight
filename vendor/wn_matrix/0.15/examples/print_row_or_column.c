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
//       Example to demonstrate how to use wn_matrix routines to read in
//       a matrix from an ascii file, then extract a row or column of the
//       matrix, output the non-zero elements to the screen, and free the
//       matrix and row or column.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include "WnMatrix.h"

int main( int argc, char *argv[] ) {

  size_t i, i_line, *a_indices;
  double *a_elements;
  WnMatrix *p_matrix;
  WnMatrix__Line *p_line;
  
  /*============================================================================
  // Check arguments.
  //==========================================================================*/

  if ( argc != 4 ) {
    fprintf(
      stderr, "\nUsage: %s in_file line i_line\n\n", argv[0]
    );
    fprintf(
      stderr, "  in_file = input matrix xmlfile\n\n"
    );
    fprintf(
      stderr, "  line = row or column\n\n"
    );
    fprintf(
      stderr, "  i_line = row or column number\n\n"
    );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Get matrix.
  //==========================================================================*/

  p_matrix = WnMatrix__new_from_xml( argv[1], NULL );

  /*============================================================================
  // Check input string.
  //==========================================================================*/

  if( strcmp( "row", argv[2] ) && strcmp( "column", argv[2] ) ) {
      fprintf( stderr, "%s should be row or column\n", argv[2] );
      return EXIT_FAILURE;
  }

  /*============================================================================
  // Get row or column.
  //==========================================================================*/

  i_line = (size_t) atoi( argv[3] );

  if( !strcmp( "row", argv[2] ) ) {
     p_line = WnMatrix__getRow( p_matrix, i_line );
     printf(
       "\nThe non-zero elements of row %lu are:\n\n",
       (unsigned long) i_line
     );
     printf( "Column        Value\n" );
     printf( "------      -------\n" );
  } else {
     p_line = WnMatrix__getColumn( p_matrix, i_line );
     printf(
       "\nThe non-zero elements of column %lu are:\n\n",
       (unsigned long) i_line
     );
     printf( "   Row        Value\n" );
     printf( "------      -------\n" );
  }

  /*============================================================================
  // Print out elements.
  //==========================================================================*/

  a_indices = WnMatrix__Line__getNonZeroIndices( p_line );
  a_elements = WnMatrix__Line__getNonZeroElements( p_line );

  for ( i = 0; i < WnMatrix__Line__getNumberOfElements( p_line ); i++ ) {

    printf( "%6lu   %10.4f\n", (unsigned long) a_indices[i], a_elements[i] );

  }

  printf( "\n" );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  WnMatrix__Line__free( p_line );
  WnMatrix__free( p_matrix );

  return EXIT_SUCCESS;

}

