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
//       store elements in a matrix, then write the elements with absolute
//       magnitude larger than the cutoff parameter to a file,
//       and finally to free the matrix.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include "WnMatrix.h"

int main( int argc, char *argv[] ) {

  double d_cutoff;
  WnMatrix *p_my_matrix;
  
  /*============================================================================
  // Check for correct arguments.
  //==========================================================================*/

  if ( argc != 3 ) {
    fprintf(
      stderr, "\nUsage: %s cutoff filename\n\n", argv[0]
    );
    fprintf(
      stderr, "  cutoff = cutoff parameter value\n\n"
    );
    fprintf(
      stderr, "  filename = data output file\n\n"
    );
    return EXIT_FAILURE;
  }

  d_cutoff = atof( argv[1] );

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
  // Write matrix to output file.
  //==========================================================================*/

  if(
      !WnMatrix__writeMatrixToAsciiFile( p_my_matrix, argv[2], d_cutoff )
  ) {
       fprintf( stderr, "Couldn't write out matrix!\n" );
  }

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  WnMatrix__free( p_my_matrix );

  return EXIT_SUCCESS;

}

