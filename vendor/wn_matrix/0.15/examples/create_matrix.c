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
//       Example to demonstrate how to use wn_matrix routines to create and
//       store elements in a matrix, print out the matrix to screen,
//       scale the matrix elements, print them out, and free the matrix.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include "WnMatrix.h"

int main( void ) {

  size_t i_row, i_col;
  WnMatrix * p_my_matrix;
  
  /*============================================================================
  // Create a 3 x 3 matrix called my_matrix and a pointer to it
  // called p_my_matrix.
  //==========================================================================*/

  p_my_matrix = WnMatrix__new( 3L, 3L );

  /* Assign the following matrix:

     | 10.  0.   3. |
     |  0.  0.   0. |
     | -5.  2.   0. |
  */

  WnMatrix__assignElement(
    p_my_matrix, 1L, 2L, 1.
  );  /* Remove this below */

  WnMatrix__assignElement( p_my_matrix, 1L, 1L, 10. );

  WnMatrix__assignElement( p_my_matrix, 3L, 1L, -2.5 );

  WnMatrix__assignElement( p_my_matrix, 3L, 2L, 200. );
  /* Update this below */

  WnMatrix__assignElement( p_my_matrix, 3L, 1L, -2.5 );

  WnMatrix__assignElement( p_my_matrix, 1L, 3L, 3. );

  WnMatrix__updateElement( p_my_matrix, 3L, 2L, 2. );

  if( WnMatrix__removeElement( p_my_matrix, 1L, 2L ) == -1 ) {
    fprintf( stderr, "Couldn't remove element!\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Print out the matrix.
  //==========================================================================*/

  printf( "\nThe elements of the matrix are:\n\n" );
  printf( "Row   Column   Value\n" );
  printf( "---   ------   -----\n" );

  for ( i_row = 1; i_row <= 3; i_row++ ) {

    for ( i_col = 1; i_col <= 3; i_col++ ) {

      fprintf( 
        stdout,
        "%3lu   %6lu   %5.1f\n", 
        (unsigned long) i_row,
        (unsigned long) i_col,
        WnMatrix__getElement( p_my_matrix, i_row, i_col ) 
      );

    }

  } 
 
  fprintf(
     stdout,
     "\nNumber of non-zero elements = %lu\n",
     (unsigned long) WnMatrix__getNumberOfElements( p_my_matrix )
  );

  /*============================================================================
  // Double the elements of the matrix and print out.
  //==========================================================================*/

  WnMatrix__scaleMatrix( p_my_matrix, 2. );

  printf( "\nThe elements of the matrix are (when doubled):\n\n" );
  printf( "Row   Column   Value\n" );
  printf( "---   ------   -----\n" );

  for ( i_row = 1; i_row <= 3; i_row++ ) {

    for ( i_col = 1; i_col <= 3; i_col++ ) {

      fprintf(
        stdout,
        "%3lu   %6lu   %5.1f\n", 
        (unsigned long) i_row,
        (unsigned long) i_col,
        WnMatrix__getElement( p_my_matrix, i_row, i_col ) 
      );

    }

  } 
 
  fprintf(
     stdout,
     "\nNumber of non-zero elements = %lu\n",
     (unsigned long) WnMatrix__getNumberOfElements( p_my_matrix )
  );

  /*============================================================================
  // Clear matrix.  Since we only clear p_my_matrix, we can still use it.
  //==========================================================================*/

  printf( "\nNow clear matrix:\n\n" );

  WnMatrix__clear( p_my_matrix );

  fprintf(
     stdout,
     "Number of non-zero elements = %lu\n\n",
     (unsigned long) WnMatrix__getNumberOfElements( p_my_matrix )
  );

  /*============================================================================
  // Assign element and print number of elements.
  //==========================================================================*/

  fprintf( stdout, "Now add element:\n\n" );

  i_row = 1; i_col = 1;
  WnMatrix__assignElement( p_my_matrix, i_row, i_col, 1. );

  fprintf(
     stdout,
     "Number of non-zero elements = %lu\n\n",
     (unsigned long) WnMatrix__getNumberOfElements( p_my_matrix )
  );

  /*============================================================================
  // Now free matrix.
  //==========================================================================*/

  fprintf( stdout, "Now free matrix.\n\n" );
  WnMatrix__free( p_my_matrix );

  /*============================================================================
  // Since p_my_matrix freed, must reallocate to reuse.
  //==========================================================================*/

  printf( "Now re-create matrix.\n\n" );
  p_my_matrix = WnMatrix__new( 3L, 3L );

  fprintf(
     stdout,
     "Number of non-zero elements = %lu\n\n",
     (unsigned long) WnMatrix__getNumberOfElements( p_my_matrix )
  );

  /*============================================================================
  // Assign element and print number of elements.
  //==========================================================================*/

  fprintf( stdout, "Now add element:\n\n" );

  i_row = 1; i_col = 1;
  WnMatrix__assignElement( p_my_matrix, i_row, i_col, 1. );

  printf(
     "Number of non-zero elements = %lu\n\n",
     (unsigned long) WnMatrix__getNumberOfElements( p_my_matrix )
  );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  WnMatrix__free( p_my_matrix );

  return EXIT_SUCCESS;

}

