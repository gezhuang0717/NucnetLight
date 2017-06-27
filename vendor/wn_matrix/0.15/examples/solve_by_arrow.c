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
//       convert the matrix to arrow form, solve the matrix equation
//       using the wn_matrix arrow solver, and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include "WnMatrix.h"

void solution_print( WnMatrix *, gsl_vector *, gsl_vector * );

int main( int argc, char **argv ) {

  WnMatrix * p_my_matrix;
  WnMatrix__Arrow *p_arrow;
  gsl_vector *p_rhs_vector, *p_solution_vector, *p_b;

  if ( argc != 3 ) {
    fprintf(
      stderr, "\nUsage: %s matrix_file vector_file \n\n", argv[0]
    );
    fprintf(
      stderr, "  matrix_file = input matrix xml data file\n\n"
    );
    fprintf(
      stderr, "  vector_file = rhs vector xml data file\n\n"
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
  // Check that rhs vector xml is valid.
  //==========================================================================*/

  if( !WnMatrix__is_valid_vector_input_xml( argv[2] ) ) {
    fprintf( stderr, "Not a valid input vector xml file!\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Get the rhs vector.
  //==========================================================================*/

  p_rhs_vector = WnMatrix__new_gsl_vector_from_xml( argv[2], NULL );

  /*============================================================================
  // Get a working copy of the rhs vector.  We do this because the arrow
  // solver modifies the right-hand-side vector and we desire a copy of the
  // original vector for later use.
  //==========================================================================*/

  p_b = gsl_vector_alloc( WnMatrix__get_gsl_vector_size( p_rhs_vector ) );
  gsl_vector_memcpy( p_b, p_rhs_vector );

  /*============================================================================
  // Convert to arrow form.  Use a wing width of 3.
  //==========================================================================*/

  p_arrow = WnMatrix__getArrow( p_my_matrix, 3L );

  /*============================================================================
  // Solve the matrix equation.
  //==========================================================================*/

  p_solution_vector = WnMatrix__Arrow__solve( p_arrow, p_b );

  /*============================================================================
  // Done with the arrow matrix and rhs copy, so free them.
  //==========================================================================*/

  gsl_vector_free( p_b );
  WnMatrix__Arrow__free( p_arrow );

  /*============================================================================
  // Print output.
  //==========================================================================*/

  solution_print( p_my_matrix, p_rhs_vector, p_solution_vector );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  gsl_vector_free( p_rhs_vector );
  gsl_vector_free( p_solution_vector );

  WnMatrix__free( p_my_matrix );

  return EXIT_SUCCESS;

}
/*##############################################################################
// solution_print().
//############################################################################*/

void
solution_print(
  WnMatrix *p_matrix,
  gsl_vector *p_rhs_vector,
  gsl_vector *p_solution_vector
)
{

  size_t i_rows;
  gsl_vector *p_check_vector;

  fprintf( stdout, "\nSolution vector:\n\n" );
  for( i_rows = 0; i_rows < WnMatrix__getNumberOfRows( p_matrix ); i_rows++ )
    fprintf(
      stdout,
      "i = %lu  a_sol[i] = %15.10e\n",
      (unsigned long) i_rows,
      gsl_vector_get( p_solution_vector, i_rows )
    );

  /*============================================================================
  // Compute matrix times solution vector.
  //==========================================================================*/

  p_check_vector =
    WnMatrix__computeMatrixTimesVector(
      p_matrix, p_solution_vector
    );

  /*============================================================================
  // Print out matrix times solution vector and compare to rhs vector.
  //==========================================================================*/

  fprintf( stdout, "\nMatrix times solution vector vs. rhs vector:\n\n" );

  for( i_rows = 0; i_rows < WnMatrix__getNumberOfRows( p_matrix ); i_rows++ )
  {
    fprintf(
      stdout,
      "i = %lu  matrix * a_sol[i] = %e  rhs[i] = %e\n",
      (unsigned long) i_rows,
      gsl_vector_get( p_check_vector, i_rows ),
      gsl_vector_get( p_rhs_vector, i_rows )
    );
  }

  fprintf( stdout, "\n" );

  gsl_vector_free( p_check_vector );

}
