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
//       a matrix from XML input, convert to arrow matrix form,
//       print out data about the arrow matrix, and free the matrices.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include "WnMatrix.h"

int main( int argc, char *argv[] ) {

  WnMatrix *p_my_matrix;
  WnMatrix__Arrow *p_arrow;
  size_t i_rows, i_band, i_wing;

  /*============================================================================
  // Check arguments.
  //==========================================================================*/

  if ( argc != 3 ) {
    fprintf(
      stderr, "\nUsage: %s filename wing_width\n\n", argv[0]
    );
    fprintf(
      stderr, "  filename = input xml data file\n\n"
    );
    fprintf(
      stderr, "  wing_width = width of arrow wings\n\n"
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

  p_my_matrix = WnMatrix__new_from_xml( argv[1], NULL );

  /*============================================================================
  // Get arrow matrix.
  //==========================================================================*/

  p_arrow = WnMatrix__getArrow( p_my_matrix, (size_t) atol( argv[2] ) );

  /*============================================================================
  // Free original matrix.
  //==========================================================================*/

  WnMatrix__free( p_my_matrix );

  /*============================================================================
  // Print arrow matrix data.
  //==========================================================================*/

  i_rows = WnMatrix__Arrow__getNumberOfRows( p_arrow );
  i_band = WnMatrix__Arrow__getBandWidth( p_arrow );
  i_wing = WnMatrix__Arrow__getWingWidth( p_arrow );

  fprintf( stdout, "\nNumber of rows = %lu\n", (unsigned long) i_rows );
  fprintf( stdout, "Wing width = %lu\n", (unsigned long) i_wing );
  fprintf( stdout, "Band width = %lu\n\n", (unsigned long) i_band );

  fprintf(
    stdout,
    "Memory for arrow matrix / memory for square matrix = %g\n\n",
    (double) (
      ( i_band + 2 * i_wing ) * ( i_rows - i_wing ) + ( i_wing * i_wing )
    ) /
    (double) ( i_rows * i_rows )
  );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  WnMatrix__Arrow__free( p_arrow );

  return EXIT_SUCCESS;

}
