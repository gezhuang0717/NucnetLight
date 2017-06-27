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
//       solve the linear matrix equation dY/dt = AY + P, where P is
//       a constant vector, with Sparskit.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <WnSparseSolve.h>
#include "solution_check.c"

int
main( int argc, char *argv[] )
{

  WnSparseSolve__Phi *p_my_phi;
  WnMatrix *p_matrix; 
  size_t i;
  double d_t;
  gsl_vector *p_initial_vector, *p_constant_vector, *p_my_solution_vector;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if( argc != 8 && argc != 9 ) {
    fprintf(
      stderr,
      "\nUsage: %s matrix_file initial_file constant_vector_file time itmax workspace tol\n\n",
      argv[0]
    );
    fprintf( stderr, "  matrix_file = name of input matrix file\n\n" );
    fprintf(
      stderr,
      "  initial_file = name of input initial vector file\n\n"
    );
    fprintf( stderr,
      "  constant_vector_file = name of input constant vector file\n\n"
    );
    fprintf( stderr, "  time = time over which to evolve\n\n" );
    fprintf( stderr, "  itmax = maximum number of iterations\n\n" );
    fprintf(
      stderr,
      "  workspace = integer size of workspace (0 < workspace < 50)\n\n"
    );
    fprintf( stderr, "  tol = tolerance (> 0)\n\n" );
    fprintf( stderr, "  debug = debug (optional)\n\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Assign time.
  //==========================================================================*/

  d_t = atof( argv[4] );

  if( d_t < 0 ) {
    printf( "Time should be > 0.\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Get matrix.
  //==========================================================================*/

  p_matrix = WnMatrix__new_from_xml( argv[1], NULL );

  /*============================================================================
  // Get initial vector.
  //==========================================================================*/

  p_initial_vector =
    WnMatrix__new_gsl_vector_from_xml( argv[2], NULL );

  /*============================================================================
  // Check that initial vector has the appropriate length.
  //==========================================================================*/

  if(
    WnMatrix__get_gsl_vector_size( p_initial_vector )
    != WnMatrix__getNumberOfColumns( p_matrix )
  ) {
    printf(
      "\nInput initial vector does not have the right number of entries!\n"
    );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Get the constant vector file.
  //==========================================================================*/

  p_constant_vector =
    WnMatrix__new_gsl_vector_from_xml( argv[3], NULL );

  /*============================================================================
  // Check that constant vector has the appropriate length.
  //==========================================================================*/

  if(
    WnMatrix__get_gsl_vector_size( p_constant_vector )
    != WnMatrix__getNumberOfColumns( p_matrix )
  ) {
    printf(
      "\nInput initial vector does not have the right number of entries!\n"
    );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Get solver.
  //==========================================================================*/

  p_my_phi = WnSparseSolve__Phi__new( );

  /*============================================================================
  // Change settings.
  //==========================================================================*/

  WnSparseSolve__Phi__updateMaximumIterations( p_my_phi, atoi( argv[5] ) );

  WnSparseSolve__Phi__updateWorkSpace( p_my_phi, atoi( argv[6] ) );

  WnSparseSolve__Phi__updateTolerance( p_my_phi, atof( argv[7] ) );

  /*============================================================================
  // Set debug.
  //==========================================================================*/

  if( argv[8] ) {
    if( strcmp( argv[8], "debug" ) == 0 )
       WnSparseSolve__Phi__setDebug( p_my_phi );
  }

  /*============================================================================
  // Get solution.
  //==========================================================================*/

  p_my_solution_vector =
    WnSparseSolve__Phi__solve(
      p_my_phi, p_matrix, p_initial_vector, p_constant_vector, d_t
    );

  /*============================================================================
  // Check that the solution succeeded.
  //==========================================================================*/

  if( !is_valid_solution( p_my_solution_vector, argv[8] ) ) {
    gsl_vector_free( p_initial_vector );
    gsl_vector_free( p_constant_vector );
    WnMatrix__free( p_matrix );
    WnSparseSolve__Phi__free( p_my_phi );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Print out solution vector.
  //==========================================================================*/

  printf( "\nSolution vector:\n\n" );
  for( i = 0; i < WnMatrix__getNumberOfRows( p_matrix ); i++ ) {
    fprintf(
      stdout,
      "i = %lu  a_w[i] = %f\n",
      (unsigned long) i,
      gsl_vector_get( p_my_solution_vector, i )
    );
  }

  /*============================================================================
  // Clean up.
  //==========================================================================*/
 
  gsl_vector_free( p_initial_vector );
  gsl_vector_free( p_constant_vector );
  gsl_vector_free( p_my_solution_vector );

  WnMatrix__free( p_matrix );

  WnSparseSolve__Phi__free( p_my_phi );

  /*============================================================================
  // Done.
  //==========================================================================*/

  return EXIT_SUCCESS;

}
