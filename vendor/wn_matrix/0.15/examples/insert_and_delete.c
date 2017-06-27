/*//////////////////////////////////////////////////////////////////////////////
// <file type="public">
//   <license>
//     This file was originally written by Joseph P. Johnson.
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
//       Example to demonstrate the use of wn_matrix routines to read in
//       an xml matrix, print out the matrix, remove a row from the matrix,
//       print the result, add a column to the matrix, and print the result.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include "WnMatrix.h"

void print_matrix( WnMatrix * );

/*##############################################################################
// main()
//############################################################################*/

int main( int argc, char *argv[] ) {

  size_t i;
  WnMatrix *p_matrix;
  gsl_vector *p_vector;

  /*============================================================================
  // Check arguments.
  //==========================================================================*/

  if ( argc != 2 ) {
    fprintf(
      stderr, "\nUsage: %s filename \n\n", argv[0]
    );
    fprintf(
      stderr, "  filename = xml matrix input file\n\n"
    );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Read xml data file.
  //==========================================================================*/

  p_matrix = WnMatrix__new_from_xml( argv[1], NULL );

  /*============================================================================
  // Print matrix.
  //==========================================================================*/

  fprintf( stdout, "Original matrix:\n" );
  print_matrix( p_matrix );

  /*============================================================================
  // Remove second row from matrix.
  //==========================================================================*/

  WnMatrix__removeRow( p_matrix, 2L );

  /*============================================================================
  // Print matrix.
  //==========================================================================*/

  fprintf( stdout, "Matrix with row removed:\n" );
  print_matrix( p_matrix );

  /*============================================================================
  // Set up gsl_vector to insert.
  //==========================================================================*/

  p_vector = gsl_vector_calloc( WnMatrix__getNumberOfRows( p_matrix ) );
  for( i = 0; i < WnMatrix__getNumberOfRows( p_matrix ); i++ ) {
    gsl_vector_set( p_vector, i, (double)( i + 1 ) );
  }

  /*============================================================================
  // Insert gsl_vector as first column of matrix.
  //==========================================================================*/

  WnMatrix__insertColumn( p_matrix, 1L, p_vector );

  /*============================================================================
  // Print matrix.
  //==========================================================================*/

  fprintf( stdout, "Matrix with column inserted:\n" );
  print_matrix( p_matrix );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  WnMatrix__free( p_matrix );
  gsl_vector_free( p_vector );

  return EXIT_SUCCESS;

}

/*##############################################################################
// print_matrix()
//############################################################################*/

void print_matrix( WnMatrix *p_matrix ) {

  size_t i_row, i_col;

  fprintf( stdout, "Row   Column   Value\n" );
  fprintf( stdout, "---   ------   -----\n" );

  for( i_row = 1; i_row <= WnMatrix__getNumberOfRows( p_matrix ); i_row++ )
  {

    for( i_col = 1; i_col <= WnMatrix__getNumberOfColumns( p_matrix ); i_col++ )
    {

      fprintf(
        stdout,
        "%3lu   %6lu   %5.1f\n",
        (unsigned long) i_row,
        (unsigned long) i_col,
        WnMatrix__getElement( p_matrix, i_row, i_col )
      );

    }

  }

  fprintf( stdout, "\n" );

  return;
}

