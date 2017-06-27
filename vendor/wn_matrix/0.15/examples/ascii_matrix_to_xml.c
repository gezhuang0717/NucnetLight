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
//       a matrix from an ascii file, store the elements, output the
//       data to an XML file, and free the matrix.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include "WnMatrix.h"

int main( int argc, char *argv[] ) {

  size_t i_rows, i_cols, i_row, i_col;
  double d_val;
  WnMatrix *p_matrix;
  WnMatrix__Coo *p_coo;
  FILE *p_in;
  
  /*============================================================================
  // Check arguments.
  //==========================================================================*/

  if ( argc != 3 && argc != 4 ) {
    fprintf(
      stderr, "\nUsage: %s in_file out_file\n\n", argv[0]
    );
    fprintf(
      stderr, "  in_file = data input file\n\n"
    );
    fprintf(
      stderr, "  out_file = output xml file\n\n"
    );
    fprintf(
      stderr,
      "  format = format code for output of value [optional--default = %%g]\n\n"
    );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Open input file.
  //==========================================================================*/

  if( ( p_in = fopen( argv[1], "r" ) ) == NULL ) {
    printf("\nCannot open input matrix file!\n");
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Determine the number of rows and columns.
  //==========================================================================*/

  i_rows = 0;
  i_cols = 0;
  while( !feof( p_in ) ) {
    fscanf(
      p_in,
      "%lu %lu %lf\n",
      (unsigned long *) &i_row,
      (unsigned long *) &i_col,
      &d_val
    );
    if( i_row > i_rows ){
      i_rows = i_row;
    }
    if( i_col > i_cols ){
      i_cols = i_col;
    }
  }

  fclose( p_in );

  /*============================================================================
  // Create a matrix.
  //==========================================================================*/

  p_matrix = WnMatrix__new( i_rows, i_cols );

  /*============================================================================
  // Reopen input file and store matrix.
  //==========================================================================*/

  p_in = fopen( argv[1], "r" );
  while( !feof( p_in ) ) {
    fscanf(
      p_in,
      "%lu %lu %lf\n",
      (unsigned long *) &i_row,
      (unsigned long *) &i_col,
      &d_val
    );
    WnMatrix__assignElement( p_matrix, (size_t) i_row, (size_t) i_col, d_val );
  }

  /*============================================================================
  // Close matrix file.
  //==========================================================================*/

  fclose( p_in );

  /*============================================================================
  // Get the coordinate matrix.
  //==========================================================================*/

  p_coo = WnMatrix__getCoo( p_matrix );

  /*============================================================================
  // Free p_matrix.
  //==========================================================================*/

  WnMatrix__free( p_matrix );

  /*============================================================================
  // Output matrix as XML.
  //==========================================================================*/

  if( argc == 3 )
    WnMatrix__Coo__writeToXmlFile( p_coo, argv[2], NULL );
  else
    WnMatrix__Coo__writeToXmlFile( p_coo, argv[2], argv[3] );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  WnMatrix__Coo__free( p_coo );

  return EXIT_SUCCESS;

}

