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
//       and store elements in a matrix and in a vector,
//       multiply the matrix by the vector,
//       output the results, multiply the transpose of the matrix by the
//       vector, output the results, and free the matrix, vector, and
//       allocated output array.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include "WnMatrix.h"

int main( int argc, char **argv ) {

  size_t i_row, i_col;
  WnMatrix * p_my_matrix;
  gsl_vector *p_input_vector, *p_output_vector;

  if ( argc != 3 ) {
    fprintf(
      stderr, "\nUsage: %s matrix_file vector_file \n\n", argv[0]
    );
    fprintf(
      stderr, "  matrix_file = input matrix xml data file\n\n"
    );
    fprintf(
      stderr, "  vector_file = input vector xml data file\n\n"
    );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Check that input matrix xml is valid.
  //==========================================================================*/

  if( !WnMatrix__is_valid_input_xml( argv[1] ) ) {
    fprintf( stderr, "Not a valid input matrix xml file!\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Get the matrix.
  //==========================================================================*/

  p_my_matrix = WnMatrix__new_from_xml( argv[1], NULL );

  /*============================================================================
  // Check that input vector xml is valid.
  //==========================================================================*/

  if( !WnMatrix__is_valid_vector_input_xml( argv[2] ) ) {
    fprintf( stderr, "Not a valid input vector xml file!\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Get the input vector.
  //==========================================================================*/

  p_input_vector = WnMatrix__new_gsl_vector_from_xml( argv[2], NULL );

  /*============================================================================
  // Multiply the matrix and vector.
  //==========================================================================*/

  printf( "Matrix times vector:\n\n" );

  p_output_vector =
    WnMatrix__computeMatrixTimesVector(
      p_my_matrix,
      p_input_vector
    );

  /*============================================================================
  // Print output.
  //==========================================================================*/

  printf( "\nRow   Column   Value   Input Vector   Output Vector\n" );
  printf( "---   ------   -----   ------------   -------------\n" );

  for(
    i_row = 1;
    i_row <= WnMatrix__getNumberOfRows( p_my_matrix );
    i_row++
  ) {

    for(
      i_col = 1;
      i_col <= WnMatrix__getNumberOfColumns( p_my_matrix );
      i_col++
    ) {

      if ( i_col == 1 ) {

        fprintf(
          stdout,
          "%3lu   %6lu   %5.1f   %12.1f   %13.1f\n",
          (unsigned long) i_row,
          (unsigned long) i_col,
          WnMatrix__getElement( p_my_matrix, i_row, i_col ),
          gsl_vector_get( p_input_vector, i_row - 1 ),
          gsl_vector_get( p_output_vector, i_row - 1 )
        );

      }
      else {

        fprintf(
          stdout,
          "%3lu   %6lu   %5.1f\n",
          (unsigned long) i_row,
          (unsigned long) i_col,
          WnMatrix__getElement( p_my_matrix, i_row, i_col )
        );
 
      }
  
    }

  }

  printf( "\n" );

  /*============================================================================
  // Done with the output vector so free it.
  //==========================================================================*/

  gsl_vector_free( p_output_vector );

  /*============================================================================
  // Multiply transpose matrix and vector.
  //==========================================================================*/

  printf( "Transpose matrix times vector:\n\n" );

  p_output_vector =
    WnMatrix__computeTransposeMatrixTimesVector(
      p_my_matrix,
      p_input_vector
    );

  /*============================================================================
  // Print output.
  //==========================================================================*/

  printf( "\nRow   Column   Value   Input Vector   Output Vector\n" );
  printf( "---   ------   -----   ------------   -------------\n" );

  for(
    i_row = 1;
    i_row <= WnMatrix__getNumberOfRows( p_my_matrix );
    i_row++
  ) {

    for(
      i_col = 1;
      i_col <= WnMatrix__getNumberOfColumns( p_my_matrix );
      i_col++
    ) {

      if ( i_col == 1 ) {

        fprintf(
          stdout,
          "%3lu   %6lu   %5.1f   %12.1f   %13.1f\n",
          (unsigned long) i_row,
          (unsigned long) i_col,
          WnMatrix__getElement( p_my_matrix, i_col, i_row ),
          gsl_vector_get( p_input_vector, i_row - 1 ),
          gsl_vector_get( p_output_vector, i_row - 1 )
        );

      }
      else {

        fprintf(
          stdout,
          "%3lu   %6lu   %5.1f\n",
          (unsigned long) i_row,
          (unsigned long) i_col,
          WnMatrix__getElement( p_my_matrix, i_col, i_row )
        );
 
      }
  
    }

  }

  printf( "\n" );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  gsl_vector_free( p_input_vector );
  gsl_vector_free( p_output_vector );
  WnMatrix__free( p_my_matrix );

  return EXIT_SUCCESS;

}
