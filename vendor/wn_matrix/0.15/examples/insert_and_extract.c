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
//       and store elements in a 2x2 matrix, insert the matrix twice into
//       a 4x4 matrix, store the 4x4 matrix in dense format, output the
//       results to a file, free the 2x2 matrix, extract a new 2x2 matrix
//       from the 4x4 matrix, output it to a file, and free the matrices.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include "WnMatrix.h"

int main( int argc, char *argv[] ) {

  size_t i, i_row, i_col;
  gsl_matrix *p_mat;
  FILE *p_outfile;
  WnMatrix *p_small_matrix, *p_large_matrix;

  /*============================================================================
  // Check arguments.
  //==========================================================================*/

  if ( argc != 2 ) {
    fprintf(
      stderr, "\nUsage: %s filename \n\n", argv[0]
    );
    fprintf(
      stderr, "  filename = data output file\n\n"
    );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Open output file.
  //==========================================================================*/

  p_outfile = fopen( argv[1], "w" );

  /*============================================================================
  // Initialize two matrices, p_small_matrix and p_large_matrix.
  //==========================================================================*/

  p_small_matrix = WnMatrix__new( 2L, 2L );

  p_large_matrix = WnMatrix__new( 4L, 4L );

  /*============================================================================
  // Assign the small matrix.
  //==========================================================================*/

  /* Assign the following matrix to p_small_matrix:

     | 1.  2. |
     | 3.  4. |
     
  */

  WnMatrix__assignElement( p_small_matrix, 1L, 1L, 1. );

  WnMatrix__assignElement( p_small_matrix, 1L, 2L, 2. );

  WnMatrix__assignElement( p_small_matrix, 2L, 1L, 3. );

  WnMatrix__assignElement( p_small_matrix, 2L, 2L, 4. );

  /*============================================================================
  // Insert the small matrix into the large one twice.
  //==========================================================================*/

  i_row = 1; i_col = 1;
  WnMatrix__insertMatrix( p_large_matrix, p_small_matrix, 1L, 1L );

  WnMatrix__insertMatrix( p_large_matrix, p_small_matrix, 3L, 3L );

  /*============================================================================
  // Store the large matrix in dense format.
  //==========================================================================*/

  p_mat = WnMatrix__getGslMatrix( p_large_matrix ); 

  /*============================================================================
  // Print out matrices.
  //==========================================================================*/

  fprintf( p_outfile, "The elements of the small matrix are:\n\n" );
  fprintf( p_outfile, "Row   Column   Value\n" );
  fprintf( p_outfile, "---   ------   -----\n" );

  for ( i_row = 1; i_row <= 2; i_row++ ) {

    for ( i_col = 1; i_col <= 2; i_col++ ) {

      fprintf(
        p_outfile,
        "%3lu   %6lu   %5.1f\n",
        (unsigned long) i_row,
        (unsigned long) i_col,
        WnMatrix__getElement( p_small_matrix, i_row, i_col )
      );

    }

  }

  fprintf( p_outfile, "\n" );

  fprintf( p_outfile, "\nThe elements of the large matrix are:\n\n" );
  fprintf( p_outfile, "Row   Column   Value\n" );
  fprintf( p_outfile, "---   ------   -----\n" );

  for ( i_row = 1; i_row <= 4; i_row++ ) {

    for ( i_col = 1; i_col <= 4; i_col++ ) {

      fprintf(
        p_outfile,
        "%3lu   %6lu   %5.1f\n",
        (unsigned long) i_row,
        (unsigned long) i_col,
        gsl_matrix_get( p_mat, i_row - 1, i_col - 1 )
      );

    }

  }

  fprintf( p_outfile, "\n" );

  /*============================================================================
  // Free the small matrix.
  //==========================================================================*/

  WnMatrix__free( p_small_matrix );

  /*============================================================================
  // Extract 2 x 2 matrix starting at row 2, column 3.
  //==========================================================================*/

  i_row = 2; i_col = 3; i = 2;
  p_small_matrix =
    WnMatrix__extractMatrix( p_large_matrix, i_row, i_col, i, i );

  /*============================================================================
  // Print out matrices.
  //==========================================================================*/

  fprintf( p_outfile, "The elements of the extracted small matrix are:\n\n" );
  fprintf( p_outfile, "Row   Column   Value\n" );
  fprintf( p_outfile, "---   ------   -----\n" );

  for ( i_row = 1;
        i_row <= WnMatrix__getNumberOfRows( p_small_matrix );
        i_row++
  ) {

    for ( i_col = 1;
          i_col <= WnMatrix__getNumberOfColumns( p_small_matrix );
          i_col++
    ) {

      fprintf(
        p_outfile,
        "%3lu   %6lu   %5.1f\n",
        (unsigned long) i_row,
        (unsigned long) i_col,
        WnMatrix__getElement( p_small_matrix, i_row, i_col )
      );

    }

  }

  fprintf( p_outfile, "\n" );

  /*============================================================================
  // Close output matrix.
  //==========================================================================*/

  fclose( p_outfile );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  WnMatrix__free( p_small_matrix );
  WnMatrix__free( p_large_matrix );

  gsl_matrix_free( p_mat );

  return EXIT_SUCCESS;

}

