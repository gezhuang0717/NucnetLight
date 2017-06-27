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
//       and store elements in a matrix, then store the matrix in compressed
//       sparse row format and output the results to a file and free
//       the matrix.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include "WnMatrix.h"

int main( int argc, char *argv[] ) {

  size_t i, i_row;
  WnMatrix * p_my_matrix;
  WnMatrix__Csr *p_csr;
  WnMatrix__Coo *p_coo;
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

  /*============================================================================
  // Open output file.
  //==========================================================================*/

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
  // Store the matrix in compressed sparse row (CSR) format and get arrays.
  //==========================================================================*/

  p_csr = WnMatrix__getCsr( p_my_matrix );


  /*============================================================================
  // Print matrix in CSR format.
  //==========================================================================*/

  fprintf( p_outfile, "Matrix in Compressed Sparse Row format:\n\n" );
  fprintf( p_outfile, "Row Pointers\n" );
  fprintf( p_outfile, "------------\n" );

  for ( i_row = 0; i_row < 4; i_row++ ) {

    fprintf(
      p_outfile,
      "%12lu\n",
      (unsigned long) ( WnMatrix__Csr__getRowPointerVector( p_csr ) )[i_row]
    );

  }

  fprintf( p_outfile, "\nIndex   Column   Value\n" );
  fprintf( p_outfile, "-----   ------   -----\n" );
  
  for(
    i = 0;
    i < WnMatrix__getNumberOfElements( p_my_matrix );
    i++
  ) {

    fprintf(
      p_outfile,
      "%5lu   %6lu   %5.1f\n",
      ( unsigned long) i,
      (unsigned long) ( WnMatrix__Csr__getColumnVector( p_csr ) )[i],
      ( WnMatrix__Csr__getValueVector( p_csr ) )[i]
    );

  }

  fprintf( p_outfile, "\n" );

  /*============================================================================
  // Done with CSR matrix, so free it.
  //==========================================================================*/

  WnMatrix__Csr__free( p_csr );

  /*============================================================================
  // Get coordinate matrix and arrays.
  //==========================================================================*/

  p_coo = WnMatrix__getCoo( p_my_matrix );

  /*============================================================================
  // Print matrix in coordinate format for comparison.
  //==========================================================================*/

  fprintf( p_outfile, "\nMatrix in coordinate format (for comparison):\n\n" );
  fprintf( p_outfile, "Row  Column  Value\n" );
  fprintf( p_outfile, "---  ------  -----\n" );

  for( i = 0; i < WnMatrix__getNumberOfElements( p_my_matrix ); i++ ) {

     fprintf(
       p_outfile,
       "%3lu%8lu%7.1f\n", 
       (unsigned long) ( WnMatrix__Coo__getRowVector( p_coo ) )[i],
       (unsigned long) ( WnMatrix__Coo__getColumnVector( p_coo ) )[i],
       ( WnMatrix__Coo__getValueVector( p_coo ) )[i]
     );

  }

  fclose( p_outfile );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  WnMatrix__Coo__free( p_coo );
  WnMatrix__free( p_my_matrix );

  return EXIT_SUCCESS;

}
