//////////////////////////////////////////////////////////////////////////////
// This file was originally written by Bradley S. Meyer.
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//! \file
//! \brief Code for solving the matrix for a zone.
////////////////////////////////////////////////////////////////////////////////

#include "matrix_solver.h"
#ifdef SPARSKIT2
#include "ilu_solvers.h"
#endif

namespace user
{

//##############################################################################
// solve_matrix_for_zone().
//##############################################################################

gsl_vector *
solve_matrix_for_zone(
  nnt::Zone& zone,
  WnMatrix * p_matrix,
  gsl_vector * p_rhs
)
{

  gsl_vector * p_sol;
  WnMatrix__Arrow * p_arrow;
#ifdef SPARSKIT2
  typedef boost::unordered_map<std::string, std::string> unordered_map;
  unordered_map param_map;
#endif

#ifndef SPARSKIT2
  if(
    zone.hasProperty( nnt::s_ITER_SOLVER ) ||
    zone.hasProperty( nnt::s_ITER_SOLVER_T9 )
  )
  {
    std::cerr << std::endl;
    std::cerr << "Code was not compiled against Sparskit2. Either recompile"
       << std::endl;
    std::cerr << "the code with Sparskit2 or remove the iterative solver zone"
       << std::endl;
    std::cerr << "properties." << std::endl << std::endl;
    exit( EXIT_FAILURE );
  }
#endif

  //============================================================================
  // Call optional matrix function.
  //============================================================================

  if( zone.getMatrixFunction() )
  {
    zone.getMatrixFunction()
    (
      p_matrix,
      p_rhs,
      zone.getMatrixFunctionData()
    );
  }

#ifdef SPARSKIT2

  //============================================================================
  // Try iterative solver, if set.
  //============================================================================

  if(
    zone.hasProperty( nnt::s_ITER_SOLVER ) &&
    zone.hasProperty( nnt::s_ITER_SOLVER_T9 ) &&
    boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) )
    <
    boost::lexical_cast<double>( zone.getProperty( nnt::s_ITER_SOLVER_T9 ) )
  )
  {

    //--------------------------------------------------------------------------
    // Get solver.  Set solver properties.
    //--------------------------------------------------------------------------

    param_map.insert(
      unordered_map::value_type(
        nnt::s_ITER_SOLVER,
        zone.getProperty( nnt::s_ITER_SOLVER )
      )
    );


    if( zone.hasProperty( nnt::s_ITER_SOLVER_MAX_ITERATIONS ) )
    {
      param_map.insert(
        unordered_map::value_type(
          nnt::s_ITER_SOLVER_MAX_ITERATIONS,
          zone.getProperty( nnt::s_ITER_SOLVER_MAX_ITERATIONS )
        )
      );
    }

    if( zone.hasProperty( nnt::s_ITER_SOLVER_REL_TOL ) )
    {
      param_map.insert(
        unordered_map::value_type(
          nnt::s_ITER_SOLVER_REL_TOL,
          zone.getProperty( nnt::s_ITER_SOLVER_REL_TOL )
        )
      );
    }
    else
    {
      param_map.insert(
        unordered_map::value_type(
          nnt::s_ITER_SOLVER_REL_TOL,
          "1.e-16"
        )
      );
    }

    if( zone.hasProperty( nnt::s_ITER_SOLVER_ABS_TOL ) )
    {
      param_map.insert(
        unordered_map::value_type(
          nnt::s_ITER_SOLVER_ABS_TOL,
          zone.getProperty( nnt::s_ITER_SOLVER_ABS_TOL )
        )
      );
    }
    else
    {
      param_map.insert(
        unordered_map::value_type(
          nnt::s_ITER_SOLVER_ABS_TOL,
          "0"
        )
      );
    }

    if( zone.hasProperty( nnt::s_ITER_SOLVER_CONVERGENCE_METHOD ) )
    {
      param_map.insert(
        unordered_map::value_type(
          nnt::s_ITER_SOLVER_CONVERGENCE_METHOD,
          zone.getProperty( nnt::s_ITER_SOLVER_CONVERGENCE_METHOD )
        )
      );
    }

    if( zone.hasProperty( nnt::s_ITER_SOLVER_DEBUG ) )
    {
      param_map.insert(
        unordered_map::value_type(
          nnt::s_ITER_SOLVER_DEBUG,
          zone.getProperty( nnt::s_ITER_SOLVER_DEBUG )
        )
      );
    }

    //--------------------------------------------------------------------------
    // Get solution.
    //--------------------------------------------------------------------------

    p_sol =
      solve_sparse_matrix_with_ilu_preconditioner(
        p_matrix,
        p_rhs,
        param_map
      );

    if( p_sol ) return p_sol;

  }

#endif

  //============================================================================
  // Normal solver.
  //============================================================================

  if(
    zone.hasProperty( nnt::s_SOLVER ) &&
    zone.getProperty( nnt::s_SOLVER ) == nnt::s_ARROW
  )
  {
    p_arrow =
      WnMatrix__getArrow(
        p_matrix,
        boost::lexical_cast<unsigned long>(
          zone.getProperty( nnt::s_ARROW_WIDTH )
        )
      );
    p_sol = WnMatrix__Arrow__solve( p_arrow, p_rhs );
    WnMatrix__Arrow__free( p_arrow );
  }
  else
  {
    p_sol = WnMatrix__solve( p_matrix, p_rhs );
  }

  return p_sol;

}

#ifdef SPARSKIT2

//##############################################################################
// solve_sparse_matrix_with_ilu_preconditioner().
//##############################################################################

gsl_vector *
solve_sparse_matrix_with_ilu_preconditioner(
  WnMatrix * p_matrix,
  gsl_vector * p_rhs,
  boost::unordered_map<std::string,std::string>& param_map
)
{

  gsl_vector * p_my_guess_vector, * p_sol;
  WnSparseSolve__Mat * p_my_solver;
  ilu_structure * p_ilu_structure;
  size_t i_delta;
  double d_drptol;
  boost::unordered_map<std::string, std::string>::iterator it;

  //============================================================================
  // Get the guess vector and zero it out.
  //============================================================================

  p_my_guess_vector =
    gsl_vector_calloc(
      WnMatrix__get_gsl_vector_size( p_rhs )
    );

  //============================================================================
  // Get solver.  Set solver properties.
  //============================================================================

  p_my_solver = WnSparseSolve__Mat__new();

  it = param_map.find( nnt::s_ITER_SOLVER );

  if( it != param_map.end() )
  {
    WnSparseSolve__Mat__updateSolverMethod(
      p_my_solver,
      it->second.c_str()
    );
  }
  else
  {
    std::cerr << "No sparse solver set." << std::endl;
    exit( EXIT_FAILURE );
  }

  it = param_map.find( nnt::s_ITER_SOLVER_MAX_ITERATIONS );

  if( it != param_map.end() )
  {
    WnSparseSolve__Mat__updateMaximumIterations(
      p_my_solver,
      boost::lexical_cast<int>( it->second )
    );
  }

  it = param_map.find( nnt::s_ITER_SOLVER_REL_TOL );

  if( it != param_map.end() )
  {
    WnSparseSolve__Mat__updateRelativeTolerance(
      p_my_solver,
      boost::lexical_cast<double>( it->second )
    );
  }

  it = param_map.find( nnt::s_ITER_SOLVER_ABS_TOL );

  if( it != param_map.end() )
  {
    WnSparseSolve__Mat__updateAbsoluteTolerance(
      p_my_solver,
      boost::lexical_cast<double>( it->second )
    );
  }

  it = param_map.find( nnt::s_ITER_SOLVER_CONVERGENCE_METHOD );

  if( it != param_map.end() )
  {
    WnSparseSolve__Mat__updateConvergenceMethod(
      p_my_solver,
      it->second.c_str()
    );
  }

  it = param_map.find( nnt::s_ITER_SOLVER_DEBUG );

  if( it != param_map.end() && it->second == "yes" )
  {
    WnSparseSolve__Mat__setDebug( p_my_solver );
  }

  //============================================================================
  // Set the preconditioner.
  //============================================================================

  it = param_map.find( nnt::s_ILU_DELTA );

  if( it != param_map.end() )
    i_delta = boost::lexical_cast<size_t>( it->second );
  else
    i_delta = 1;

  it = param_map.find( nnt::s_ILU_DROP_TOL );

  if( it != param_map.end() )
    d_drptol = boost::lexical_cast<double>( it->second );
  else
    d_drptol = 0.;

  WnSparseSolve__Mat__updatePreconditionerSolver(
    p_my_solver,
    (WnSparseSolve__Mat__PreconditionerSolver) ilu_solver
  );

  WnSparseSolve__Mat__updatePreconditionerTransposeSolver(
    p_my_solver,
    (WnSparseSolve__Mat__PreconditionerTransposeSolver)
      ilu_transpose_solver
  );

  p_ilu_structure =
    ilu_preconditioner_new(
      p_matrix,
      i_delta,
      d_drptol
    );

  WnSparseSolve__Mat__updatePreconditionerUserData(
    p_my_solver, p_ilu_structure
  );

  //============================================================================
  // Solve the matrix equation.
  //============================================================================

  p_sol =
    WnSparseSolve__Mat__solve(
      p_my_solver,
      p_matrix,
      p_rhs,
      p_my_guess_vector
    );

  //============================================================================
  // Clean up and return if valid solution.
  //============================================================================

  WnSparseSolve__Mat__free( p_my_solver );
  ilu_preconditioner_free( p_ilu_structure );
  gsl_vector_free( p_my_guess_vector );

  return p_sol;

}

//##############################################################################
// phi__solve__parallel().
//##############################################################################

gsl_vector *
phi__solve__parallel(
  WnSparseSolve__Phi *self,
  std::vector<WnMatrix *>& matrices,
  const gsl_vector *p_input_vector,
  const gsl_vector *p_const_vector,
  double d_tn
)
{

  size_t i, i_rows, n_fpar = 6;
  int i_indic, i_ierr, i_iter, i_workspace;
  int i_ipar[5] = {0,0,0,0,0};
  double d_tol;
  gsl_vector *p_x, *p_y, *p_u, *p_output_vector, *p_fpar;

  //============================================================================
  // Check input.
  //============================================================================

  if( matrices.size() == 0 || !p_input_vector || !p_const_vector )
    WNSPARSESOLVE__ERROR( "Invalid input" );

  //============================================================================
  // Find debugging information.
  //============================================================================

  if( self->iDebug == 1 )
      printf( "\nUsing phipro\n\n" );

  //============================================================================
  // Set debugging flag.
  //============================================================================

  i_ipar[3] = self->iDebug;

  //============================================================================
  // Get number of rows.
  //============================================================================

  i_rows = WnMatrix__getNumberOfRows( matrices[0] );

  //============================================================================
  // Set local variables from structure.
  //============================================================================

  i_workspace = self->iWorkSpace;
  d_tol = self->dTolerance;

  //============================================================================
  // Vector to handle multiplication results.
  //============================================================================

  std::vector<gsl_vector *> xx( matrices.size() );

  //============================================================================
  // Allocate memory for vectors.
  //============================================================================

  p_x = gsl_vector_calloc( i_rows );
  p_y = gsl_vector_calloc( i_rows );
  p_u = gsl_vector_calloc( i_rows * ( (size_t) i_workspace + 1 ) );
  p_output_vector = gsl_vector_alloc( p_input_vector->size );
  p_fpar = gsl_vector_calloc( n_fpar );

  if( !p_x || !p_y || !p_u || !p_output_vector || !p_fpar )
      WNSPARSESOLVE__ERROR( "Couldn't allocate memory in WnSparseSolve" );

  gsl_vector_memcpy( p_output_vector, p_input_vector );

  i_indic = 0;
  i_iter = 0;
  i_ierr = 0;

  while( i_indic != 1 ) {

    if( i_iter == self->iIterMax ) {
      if( self->iDebug == 1 )
        fprintf( stdout, "Exceeded maximum number of iterations.\n" );
      gsl_vector_free( p_x );
      gsl_vector_free( p_y );
      gsl_vector_free( p_u );
      gsl_vector_free( p_output_vector );
      gsl_vector_free( p_fpar );
      return NULL;
    }

    phipro_(
      &i_rows,
      &i_workspace,
      &d_tol,
      &d_tn,
      p_output_vector->data,
      p_const_vector->data,
      p_u->data,
      p_x->data,
      p_y->data,
      p_fpar->data,
      i_ipar,
      &i_indic,
      &i_ierr
    );

    gsl_vector_free( p_y );

    #pragma omp parallel for schedule( dynamic, 1 )
      for( i = 0; i < matrices.size(); i++ )
      {
        xx[i] =
          WnMatrix__computeMatrixTimesVector(
            matrices[i],
            p_x
          );
      }

    size_t i_add = 2;

    size_t i_number = matrices.size();

    while( i_number / 2 > 0 )
    {

      size_t i_c = i_add / 2;

      if( GSL_IS_ODD( i_number ) )
      {
        gsl_vector_add(
          xx[(i_number/2)*i_add - i_c],
          xx[(i_number/2)*i_add]
        );
        gsl_vector_free( xx[(i_number/2)*i_add] );
      }

      i_number /= 2;

      #pragma omp parallel for schedule( dynamic, 1 )
        for( i = 0; i < i_number; i++ )
        {
          gsl_vector_add( xx[i*i_add], xx[i*i_add+i_c] );
          gsl_vector_free( xx[i*i_add+i_c] );
        }

      i_add *= 2;

    }

    p_y = gsl_vector_calloc( i_rows );

    gsl_vector_memcpy( p_y, xx[0] );

    gsl_vector_free( xx[0] );

    gsl_vector_scale( p_y, -1. );

    ++i_iter;

  }

  if( self->iDebug == 1 )
    printf( "%d iterations.\n", i_iter );

  //============================================================================
  // Clean up and exit.
  //============================================================================

  gsl_vector_free( p_x );
  gsl_vector_free( p_y );
  gsl_vector_free( p_u );
  gsl_vector_free( p_fpar );

  return p_output_vector;

}

#endif  // SPARSKIT2

} //namespace user
