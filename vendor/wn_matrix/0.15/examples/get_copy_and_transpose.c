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
//       a matrix from an ascii file, get a copy of the matrix and the
//       transpose, output the elements to a file, and then free all
//       matrices and their allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include "WnMatrix.h"

int main( int argc, char *argv[] ) {

  size_t i_row, i_col;
  WnMatrix *p_matrix, *p_copy, *p_transpose;
  FILE *p_outfile;
  
  /*============================================================================
  // Check arguments.
  //==========================================================================*/

  if ( argc != 3 ) {
    fprintf(
      stderr, "\nUsage: %s in_file out_file\n\n", argv[0]
    );
    fprintf(
      stderr, "  in_file = input matrix xml file\n\n"
    );
    fprintf(
      stderr, "  out_file = data output file\n\n"
    );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Get matrix.
  //==========================================================================*/

  p_matrix = WnMatrix__new_from_xml( argv[1], NULL );

  /*============================================================================
  // Open output file.
  //==========================================================================*/

  if( ( p_outfile = fopen( argv[2], "w" ) ) == NULL ) {
    printf("\nCannot open output matrix file!\n");
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Get copy of matrix.
  //==========================================================================*/

  p_copy = WnMatrix__getCopy( p_matrix );

  /*============================================================================
  // Get transpose of matrix.
  //==========================================================================*/

  p_transpose = WnMatrix__getTranspose( p_matrix );

  /*============================================================================
  // Print output to file.
  //==========================================================================*/

  fprintf( p_outfile, "The elements of the matrices are:\n\n" );
  fprintf( p_outfile, "Row   Column    Matrix        Copy     Transpose\n" );
  fprintf( p_outfile, "---   ------    ------      -------    ---------\n" );

  for ( i_row = 1; i_row <= WnMatrix__getNumberOfRows( p_copy ); i_row++ ) {

    for ( i_col = 1; 
          i_col <= WnMatrix__getNumberOfColumns( p_copy ); 
          i_col++ 
    ) {

      fprintf(
        p_outfile,
        "%3lu   %6lu   %7.3f      %7.3f      %7.3f\n",
        (unsigned long) i_row,
        (unsigned long) i_col,
        WnMatrix__getElement( p_matrix, i_row, i_col ),
        WnMatrix__getElement( p_copy, i_row, i_col ),
        WnMatrix__getElement( p_transpose, i_row, i_col )
      );

    }

  }

  /*============================================================================
  // Close file, clean up, and exit.
  //==========================================================================*/

  fclose( p_outfile );

  WnMatrix__free( p_matrix );
  WnMatrix__free( p_copy );
  WnMatrix__free( p_transpose );

  return EXIT_SUCCESS;

}

