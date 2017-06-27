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
//     Please see the src/README.txt file in this distribution for more
//     information.
//   </license>
//   <description>
//     <abstract>
//       Example to demonstrate how to use wn_sparse_solve routines to 
//       solve the linear matrix equation dY/dt = AY with Sparskit.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <WnSparseSolve.h>
#include "print_out.c"
#include "solution_check.c"

void print_out( WnMatrix *, gsl_vector * );

int main( int argc, char *argv[] ) {

  WnSparseSolve__Exp *p_my_exp;
  WnMatrix *p_matrix; 
  gsl_vector *p_input_vector, *p_my_solution_vector;
  double d_t;
  size_t i;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if( argc != 7 && argc != 8 ) {
    fprintf(
      stderr,
      "\nUsage: %s matrix_file initial_file time itmax workspace tol\n\n",
      argv[0]
    );
    fprintf( stderr, "  matrix_file = name of input matrix file\n\n" );
    fprintf( stderr, "  initial_file = name of input initial vector file\n\n" );
    fprintf( stderr, "  time = time over which to evolve\n\n" );
    fprintf( stderr, "  itmax = maximum number of iterations\n\n" );
    fprintf( stderr,
      "  workspace: integer size of workspace (0 < workspace < 50)\n\n"
    );
    fprintf( stderr, "  tol = tolerance (> 0)\n\n" );
    fprintf( stderr, "  debug = debug (optional)\n\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Assign time, workspace, tolerance, and debug information.
  //==========================================================================*/

  d_t = atof( argv[3] );

  if( d_t < 0 ) {
    printf( "Time should be > 0.\n" );
    return 1;
  }

  /*============================================================================
  // Get the matrix.
  //==========================================================================*/

  p_matrix = WnMatrix__new_from_xml( argv[1], NULL );

  /*============================================================================
  // Get the initial vector file.
  //==========================================================================*/

  p_input_vector =
    WnMatrix__new_gsl_vector_from_xml( argv[2], NULL );

  /*============================================================================
  // Check that initial vector has the appropriate length.
  //==========================================================================*/

  if(
      WnMatrix__get_gsl_vector_size( p_input_vector )
      != WnMatrix__getNumberOfColumns( p_matrix )
  ) {
    fprintf(
      stderr,
      "\nInput initial vector does not have the right number of entries!\n"
    );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Print out.
  //==========================================================================*/

  print_out( p_matrix, p_input_vector );

  /*============================================================================
  // Get solver.
  //==========================================================================*/

  p_my_exp = WnSparseSolve__Exp__new( );

  /*============================================================================
  // Change settings.
  //==========================================================================*/

  WnSparseSolve__Exp__updateMaximumIterations( p_my_exp, atoi( argv[4] ) );

  WnSparseSolve__Exp__updateWorkSpace( p_my_exp, atoi( argv[5] ) );

  WnSparseSolve__Exp__updateTolerance( p_my_exp, atof( argv[6] ) );

  /*============================================================================
  // Set debug.
  //==========================================================================*/

  if( argv[7] ) {
    if( strcmp( argv[7], "debug" ) == 0 )
       WnSparseSolve__Exp__setDebug( p_my_exp );
  }

  /*============================================================================
  // Find solution.
  //==========================================================================*/

  p_my_solution_vector =
    WnSparseSolve__Exp__solve(
      p_my_exp,
      p_matrix,
      p_input_vector,
      d_t
    );
      
  /*============================================================================
  // Check that the solution succeeded.
  //==========================================================================*/

  if( !is_valid_solution( p_my_solution_vector, argv[7] ) ) {
    gsl_vector_free( p_input_vector );
    WnMatrix__free( p_matrix );
    WnSparseSolve__Exp__free( p_my_exp );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Print out solution vector.
  //==========================================================================*/

  printf( "\nSolution vector:\n\n" );
  for( i = 0; i < WnMatrix__getNumberOfRows( p_matrix ); i++ ) {
    fprintf(
      stdout,
      "i = %lu  sol[i] = %f\n",
      (unsigned long) i,
      gsl_vector_get( p_my_solution_vector, i )
    );
  }

  printf( "\n" );

  /*============================================================================
  // Clean up.
  //==========================================================================*/
 
  gsl_vector_free( p_input_vector );
  gsl_vector_free( p_my_solution_vector );

  WnMatrix__free( p_matrix );

  WnSparseSolve__Exp__free( p_my_exp );

  /*============================================================================
  // Done.
  //==========================================================================*/

  return EXIT_SUCCESS;

}
