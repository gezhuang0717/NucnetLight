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
//       solve a matrix equation with the matrix stored in native wn_matrix
//       format.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/
 
#include <WnSparseSolve.h>
#include "print_out.c"
#include "solution_check.c"

void print_out( WnMatrix *, gsl_vector * );

int main( int argc, char *argv[] ) {
  WnMatrix *p_matrix;
  gsl_vector *p_my_rhs_vector, *p_my_guess_vector, *p_my_solution_vector;
  WnSparseSolve__Mat *p_my_solver;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if( argc != 7 && argc != 8 ) {
    fprintf(
      stderr,
      "\nUsage: %s matrix_file rhs_file solver itmax rel_tol abs_tol debug\n\n",
      argv[0]
    );
    fprintf( stderr, "  matrix_file = name of input matrix file\n\n" );
    fprintf( stderr, "  rhs_file = name of input rhs vector file\n\n" );
    fprintf( stderr, "  solver = name of solver to use\n\n" );
    fprintf( stderr, "  itmax = maximum number of iterations\n\n" );
    fprintf( stderr, "  rel_tol = relative tolerance\n\n" );
    fprintf( stderr, "  abs_tol = absolute tolerance\n\n" );
    fprintf( stderr, "  debug = use debug mode (optional) \n\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Open matrix file.
  //==========================================================================*/

  p_matrix = WnMatrix__new_from_xml( argv[1], NULL );

  /*============================================================================
  // Get the rhs vector.
  //==========================================================================*/

  p_my_rhs_vector =
    WnMatrix__new_gsl_vector_from_xml( argv[2], NULL );

  /*============================================================================
  // Check that rhs vector has the appropriate length.
  //==========================================================================*/

  if(
      WnMatrix__get_gsl_vector_size( p_my_rhs_vector )
      != WnMatrix__getNumberOfColumns( p_matrix )
  )
  {
    fprintf(
      stderr,
      "\nInput rhs vector does not have the right number of entries!\n"
    );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Get the guess vector and zero it out.
  //==========================================================================*/

  p_my_guess_vector =
    gsl_vector_calloc( 
      WnMatrix__get_gsl_vector_size( p_my_rhs_vector )
    );

  /*============================================================================
  // Print out matrix and rhs vector.
  //==========================================================================*/

  print_out( p_matrix, p_my_rhs_vector );

  /*============================================================================
  // Get solver.
  //==========================================================================*/

  p_my_solver = WnSparseSolve__Mat__new();

  WnSparseSolve__Mat__updateSolverMethod( p_my_solver, argv[3] );

  WnSparseSolve__Mat__updateMaximumIterations( p_my_solver, atoi( argv[4] ) );

  WnSparseSolve__Mat__updateRelativeTolerance( p_my_solver, atof( argv[5] ) );

  WnSparseSolve__Mat__updateAbsoluteTolerance( p_my_solver, atof( argv[6] ) );

  if( argv[7] ) {
    if( strcmp( argv[7], "debug" ) == 0 )
      WnSparseSolve__Mat__setDebug( p_my_solver );
  }

  /*============================================================================
  // Print out which solver.
  //==========================================================================*/

  fprintf(
    stdout, 
    "\n\nUsing %s solver method.\n\n",
    WnSparseSolve__Mat__getSolverMethod( p_my_solver )
  );

  /*============================================================================
  // Solve the matrix equation.
  //==========================================================================*/

  p_my_solution_vector =
    WnSparseSolve__Mat__solve(
      p_my_solver,
      p_matrix,
      p_my_rhs_vector,
      p_my_guess_vector
    );

  /*============================================================================
  // Check that the solution succeeded.
  //==========================================================================*/

  if( !is_valid_solution( p_my_solution_vector, argv[7] ) ) {
    gsl_vector_free( p_my_rhs_vector );
    gsl_vector_free( p_my_guess_vector );
    WnMatrix__free( p_matrix );
    WnSparseSolve__Mat__free( p_my_solver );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Print solution.
  //==========================================================================*/

  solution_print( p_matrix, p_my_rhs_vector, p_my_solution_vector );

  /*============================================================================
  // Done with the solution vector.
  //==========================================================================*/

  gsl_vector_free( p_my_solution_vector );

  /*============================================================================
  // Solve the matrix equation again but with difference convergence method.
  //==========================================================================*/

  fprintf( stdout, "\nNow use rhs convergence method:\n\n" );

  WnSparseSolve__Mat__updateConvergenceMethod( p_my_solver, "rhs" );

  p_my_solution_vector =
    WnSparseSolve__Mat__solve(
      p_my_solver,
      p_matrix,
      p_my_rhs_vector,
      p_my_guess_vector
    );

  /*============================================================================
  // Check that the solution succeeded.
  //==========================================================================*/

  if( !is_valid_solution( p_my_solution_vector, argv[7] ) ) {
    gsl_vector_free( p_my_rhs_vector );
    gsl_vector_free( p_my_guess_vector );
    WnMatrix__free( p_matrix );
    WnSparseSolve__Mat__free( p_my_solver );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Print out the solution.
  //==========================================================================*/

  solution_print( p_matrix, p_my_rhs_vector, p_my_solution_vector );

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  gsl_vector_free( p_my_solution_vector );
  gsl_vector_free( p_my_rhs_vector );
  gsl_vector_free( p_my_guess_vector );

  WnMatrix__free( p_matrix );

  WnSparseSolve__Mat__free( p_my_solver );

  /*============================================================================
  // Done.
  //==========================================================================*/

  return EXIT_SUCCESS;

}
