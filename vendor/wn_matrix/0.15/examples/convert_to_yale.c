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
//       and store elements in a matrix, store the elements to Yale sparse
//       matrix format, output the results to a file, and then free the
//       matrix and arrays.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include "WnMatrix.h"

int main( int argc, char *argv[] ) {

  size_t i, i_row, i_col;
  size_t *a_ija;
  double *a_val;
  WnMatrix * p_my_matrix;
  WnMatrix__Yale *p_yale;
  FILE *p_outfile;
  
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

  p_outfile = fopen( argv[1], "w" );

  /*============================================================================
  // Create a 3 x 3 matrix.
  //==========================================================================*/

  p_my_matrix = WnMatrix__new( 3L, 3L );

  /* Assign the following matrix:

     | 10.  0.   3. |
     |  0.  0.   0. |
     | -5.  2.   0. |
  */

  WnMatrix__assignElement( p_my_matrix, 1L, 1L, 10. );

  WnMatrix__assignElement( p_my_matrix, 3L, 1L, -5. );

  WnMatrix__assignElement( p_my_matrix, 3L, 2L, 2. );

  WnMatrix__assignElement( p_my_matrix, 1L, 3L, 3. );

  /*============================================================================
  // Store matrix in Yale form and arrays.
  //==========================================================================*/

  p_yale = WnMatrix__getYale( p_my_matrix );

  a_ija = WnMatrix__Yale__getPointerVector( p_yale );
  a_val = WnMatrix__Yale__getValueVector( p_yale );

  /*============================================================================
  // Print matrix in Yale form.
  //==========================================================================*/

  fprintf( p_outfile, "Matrix in Yale Sparse format:\n\n" );

  fprintf( p_outfile, "Index   ija       sa\n" );
  fprintf( p_outfile, "-----   ---   ----------\n" );


  for ( i = 0; i < a_ija[0]; i++)
    fprintf(
      p_outfile,
      "  %lu      %lu  %10.4f\n",
      (unsigned long) i,
      (unsigned long) a_ija[i],
      a_val[i]
    ); 

  for( i = a_ija[0]; i < a_ija[a_ija[0]-1]; i++ )
    fprintf(
      p_outfile,
      "  %lu      %lu  %10.4f\n",
      (unsigned long) i,
      (unsigned long) a_ija[i],
      a_val[i]
    ); 

  /*============================================================================
  // Print matrix in dense format for comparison.
  //==========================================================================*/

  fprintf( p_outfile, "\nMatrix in dense format (for comparison):\n\n" );
  fprintf( p_outfile, "Row  Column  Value\n" );
  fprintf( p_outfile, "---  ------  -----\n" );

  for( i_row = 1; i_row <=3; i_row++ ) {

     for( i_col = 1; i_col <= 3; i_col++ ) {

         fprintf(
           p_outfile,
           "%3lu%8lu%7.1f\n", 
           (unsigned long) i_row,
           (unsigned long) i_col,
           WnMatrix__getElement( p_my_matrix, i_row, i_col )
         );
     }

  }

  fclose( p_outfile );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  WnMatrix__Yale__free( p_yale );
  WnMatrix__free( p_my_matrix );

  return EXIT_SUCCESS;

}
