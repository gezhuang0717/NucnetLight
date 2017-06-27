/*//////////////////////////////////////////////////////////////////////////////
// <file type="public">
//   <license>
//      Copyright (c) 2006-2012 Clemson University.
//
//      This file is part of the Clemson Webnucleo group's
//      wn_sparse_solve module, originally developed by Bradley S. Meyer
//      and Lih-Sin The.  For more information,
//      please see http://www.webnucleo.org.
//
//      This is free software; you can redistribute it and/or modify it
//      under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//
//      This software is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//
//      You should have received a copy of the GNU General Public License
//      along with this software (please see the "gnu_gpl.txt" file in the doc/
//      directory of this distribution); if not, write to the Free Software
//      Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
//      USA
//   </license>
//   <description>
//     <abstract>
//        Source code for the wn_sparse_solve module.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include "WnSparseSolve.h"

/*##############################################################################
// WnSparseSolve__Mat__new().
//############################################################################*/

WnSparseSolve__Mat *
WnSparseSolve__Mat__new( void )
{

  WnSparseSolve__Mat *self;

  self = (WnSparseSolve__Mat *) malloc( sizeof( WnSparseSolve__Mat ) ); 

  self->iSolverMethod = DEFAULT_SOLVER;
  self->iConvergenceMethod = DEFAULT_CONVERGE;
  self->iDebug = 0; 
  self->iIterMax = MAT_ITER_MAX;
  self->dAbsTolerance = D_ABS_TOL;
  self->dRelTolerance = D_ABS_TOL;
  self->pfPreconditionerSolve = WnSparseSolve__Mat__defaultSolve;
  self->pfPreconditionerTransposeSolve =
    WnSparseSolve__Mat__defaultTransposeSolve;
  self->pPreconditionerUserData = NULL;
  self->pfConvergenceTester = NULL;
  self->pConvergenceTesterUserData = NULL;

  return self;

}

/*##############################################################################
// WnSparseSolve__Mat__setDebug().
//############################################################################*/

void
WnSparseSolve__Mat__setDebug( WnSparseSolve__Mat *self )
{

  if( !self ) 
    WNSPARSESOLVE__ERROR( "Invalid input" );

  self->iDebug = 1;

}

/*##############################################################################
// WnSparseSolve__Mat__clearDebug().
//############################################################################*/

void
WnSparseSolve__Mat__clearDebug( WnSparseSolve__Mat *self )
{

  if( !self ) 
    WNSPARSESOLVE__ERROR( "Invalid input" );

  self->iDebug = 0;

}

/*##############################################################################
// WnSparseSolve__Mat__updateMaximumIterations().
//############################################################################*/

void
WnSparseSolve__Mat__updateMaximumIterations(
  WnSparseSolve__Mat *self, int i_iter_max
)
{

  if( !self ) 
    WNSPARSESOLVE__ERROR( "Invalid input" );

  if( i_iter_max <= 0 )
    WNSPARSESOLVE__ERROR( "Maximum number of iterations must be > 0" );

  self->iIterMax = i_iter_max;

}

/*##############################################################################
// WnSparseSolve__Mat__updateRelativeTolerance().
//############################################################################*/

void
WnSparseSolve__Mat__updateRelativeTolerance(
  WnSparseSolve__Mat *self, double d_rel_tol
)
{

  if( !self ) 
    WNSPARSESOLVE__ERROR( "Invalid input" );

  if( d_rel_tol < 0. )
    WNSPARSESOLVE__ERROR( "Relative tolerance must be non-negative" );

  self->dRelTolerance = d_rel_tol;

}

/*##############################################################################
// WnSparseSolve__Mat__updateAbsoluteTolerance().
//############################################################################*/

void
WnSparseSolve__Mat__updateAbsoluteTolerance(
  WnSparseSolve__Mat *self, double d_abs_tol
)
{

  if( !self ) 
    WNSPARSESOLVE__ERROR( "Invalid input" );

  if( d_abs_tol < 0. )
    WNSPARSESOLVE__ERROR( "Absolute tolerance must be non-negative" );

  self->dAbsTolerance = d_abs_tol;

}

/*##############################################################################
// WnSparseSolve__Mat__updatePreconditionerSolver().
//############################################################################*/

void
WnSparseSolve__Mat__updatePreconditionerSolver(
  WnSparseSolve__Mat *self,
  WnSparseSolve__Mat__PreconditionerSolver pfFunction
)
{

  if( !self ) WNSPARSESOLVE__ERROR( "Invalid input" );

  self->pfPreconditionerSolve = pfFunction;

}

/*##############################################################################
// WnSparseSolve__Mat__updatePreconditionerTransposeSolver().
//############################################################################*/

void
WnSparseSolve__Mat__updatePreconditionerTransposeSolver(
  WnSparseSolve__Mat *self,
  WnSparseSolve__Mat__PreconditionerTransposeSolver pfFunction
)
{

  if( !self ) WNSPARSESOLVE__ERROR( "Invalid input" );

  self->pfPreconditionerTransposeSolve = pfFunction;

}

/*##############################################################################
// WnSparseSolve__Mat__updateConvergenceTester().
//############################################################################*/

void
WnSparseSolve__Mat__updateConvergenceTester(
  WnSparseSolve__Mat *self,
  WnSparseSolve__Mat__ConvergenceTester pfFunction
)
{

  if( !self ) WNSPARSESOLVE__ERROR( "Invalid input" );

  self->pfConvergenceTester = pfFunction;

}

/*##############################################################################
// WnSparseSolve__Mat__updatePreconditionerUserData().
//############################################################################*/

void
WnSparseSolve__Mat__updatePreconditionerUserData(
  WnSparseSolve__Mat *self, void * p_user_data
)
{

  if( !self ) WNSPARSESOLVE__ERROR( "Invalid input" );

  self->pPreconditionerUserData = p_user_data;

}

/*##############################################################################
// WnSparseSolve__Mat__updateConvergenceTesterUserData().
//############################################################################*/

void
WnSparseSolve__Mat__updateConvergenceTesterUserData(
  WnSparseSolve__Mat *self, void * p_user_data
)
{

  if( !self ) WNSPARSESOLVE__ERROR( "Invalid input" );

  self->pConvergenceTesterUserData = p_user_data;

}

/*##############################################################################
// WnSparseSolve__Mat__updateSolverMethod().
//############################################################################*/

void
WnSparseSolve__Mat__updateSolverMethod(
  WnSparseSolve__Mat *self, const char *s_solver
)
{

  if( !strcmp( s_solver, S_CG ) ) {
    self->iSolverMethod = CG;
    return;
  }
  else if( !strcmp( s_solver, S_CGNR ) ) {
    self->iSolverMethod = CGNR;
    return;
  }
  else if( !strcmp( s_solver, S_BCG ) ) {
    self->iSolverMethod = BCG;
    return;
  }
  else if( !strcmp( s_solver, S_DBCG ) ) {
    self->iSolverMethod = DBCG;
    return;
  }
  else if( !strcmp( s_solver, S_BCGSTAB ) ) {
    self->iSolverMethod = BCGSTAB;
    return;
  }
  else if( !strcmp( s_solver, S_TFQMR ) ) {
    self->iSolverMethod = TFQMR;
    return;
  }
  else if( !strcmp( s_solver, S_FOM ) ) {
    self->iSolverMethod = FOM;
    return;
  }
  else if( !strcmp( s_solver, S_GMRES ) ) {
    self->iSolverMethod = GMRES;
    return;
  }
  else if( !strcmp( s_solver, S_FGMRES ) ) {
    self->iSolverMethod = FGMRES;
    return;
  }
  else if( !strcmp( s_solver, S_DQGMRES ) ) {
    self->iSolverMethod = DQGMRES;
    return;
  }
  else {
    WNSPARSESOLVE__ERROR( "No such solver" );
  }

}

/*##############################################################################
// WnSparseSolve__Mat__updateConvergenceMethod().
//############################################################################*/

void
WnSparseSolve__Mat__updateConvergenceMethod(
  WnSparseSolve__Mat *self, const char *s_solver
)
{

  if( !strcmp( s_solver, S_CONVERGE_INITIAL ) ) {
    self->iConvergenceMethod = CONVERGE_INITIAL ;
    return;
  }
  else if( !strcmp( s_solver, S_CONVERGE_RHS ) ) {
    self->iConvergenceMethod = CONVERGE_RHS;
    return;
  }
  else {
    WNSPARSESOLVE__ERROR( "No such convergence method" );
  }

}

/*##############################################################################
// WnSparseSolve__Mat__getSolverMethod().
//############################################################################*/

const char *
WnSparseSolve__Mat__getSolverMethod( WnSparseSolve__Mat *self )
{

  switch( self->iSolverMethod ) {

    case CG: return S_CG;
    case CGNR: return S_CGNR;
    case BCG: return S_BCG;
    case DBCG: return S_DBCG;
    case BCGSTAB: return S_BCGSTAB;
    case TFQMR: return S_TFQMR;
    case FOM: return S_FOM;
    case GMRES: return S_GMRES;
    case FGMRES: return S_FGMRES;
    case DQGMRES: return S_DQGMRES;
    default: return NULL;

  }

}

/*##############################################################################
// WnSparseSolve__Mat__solve().
//############################################################################*/

gsl_vector *
WnSparseSolve__Mat__solve(
  WnSparseSolve__Mat *self,
  WnMatrix *p_matrix,
  gsl_vector *p_rhs_vector,
  gsl_vector *p_guess_vector
) {

  size_t i, i_rows, i_elements;
  int i_iter = 0;
  int i_krylov = 16;
  int    a_ipar[16];
  double a_fpar[16], *a_wk, d_residual;
  gsl_vector *p_solution_vector, *p_in, *p_out = NULL;
  gsl_vector *p_delta, *p_working_solution;
  gsl_vector *p_tmp;

  /*============================================================================
  // Print out solver information, if desired.
  //==========================================================================*/

  if( self->iDebug == 1 ) {
    fprintf(
      stdout,
      "\n\nUsing %s solver.\n",
      WnSparseSolve__Mat__getSolverMethod( self )
    );
  }

  /*============================================================================
  // Get matrix information.
  //==========================================================================*/

  i_rows = WnMatrix__getNumberOfRows( p_matrix );
  i_elements = WnMatrix__getNumberOfElements( p_matrix );

  /*============================================================================
  // Initialize parameters.
  //==========================================================================*/

  for( i = 0; i < 16; i++ ) {
    a_ipar[i] = 0;
    a_fpar[i] = 0.;
  }

  /*============================================================================
  // Assign parameters.
  //==========================================================================*/

  a_ipar[1] = 1;
  a_ipar[4] = i_krylov;
  a_ipar[5] = self->iIterMax; 
  a_fpar[0] = self->dRelTolerance;
  a_fpar[1] = self->dAbsTolerance;

  if( self->pfConvergenceTester ) {
    switch( self->iSolverMethod )
    {
       case CG: case FOM: case GMRES: case FGMRES: case DQGMRES:
         if( self->iDebug == 1 ) {
           fprintf(
             stderr,
             "Invalid solver with user-defined convergence tester.\n"
           );
         }
         return NULL;

       default:
         break;
    }
    a_ipar[2] = 999;
  }
  else {
    switch( self->iConvergenceMethod ) {
      case 0:
        a_ipar[2] = 1;
        break;

      case 1:
        a_ipar[2] = 2;
        break;

      default:
        WNSPARSESOLVE__ERROR( "No such convergence method" );

    }
  }

  if( i_rows == 0 ) {
    WNSPARSESOLVE__ERROR( "Zero rows in matrix" );
  }

  switch( self->iSolverMethod )
  {

    case CG:
      a_ipar[3] = 5 * (int) i_rows;
      break;

    case CGNR:
      a_ipar[3] = 5 * (int) i_rows;
      break;

    case BCG:
      a_ipar[3] = 7 * (int) i_rows;
      break;

    case DBCG:
      a_ipar[3] = 11 * (int) i_rows;
      break;

    case BCGSTAB:
      a_ipar[3] =  8 * (int) i_rows;
      break;

    case TFQMR:
      a_ipar[3] = 11 * (int) i_rows;
      break;

    case FOM:
      a_ipar[3] = ( (int) i_rows + 3 ) * (i_krylov + 2 ) +
      (i_krylov + 1 ) * i_krylov / 2;
      break;

    case GMRES:
      a_ipar[3] = ( (int) i_rows + 3 ) * (i_krylov + 2 ) +
      (i_krylov + 1 ) * i_krylov / 2;
      break;

    case FGMRES:
      a_ipar[3] =   2 * (int) i_rows * (i_krylov + 1 ) +
      (i_krylov + 1 ) * i_krylov / 2 + 3 * i_krylov + 2 ;
      break;

    case DQGMRES:
      a_ipar[3] = (int) i_rows + ( a_ipar[4] + 1 ) * ( 2 * (int) i_rows + 4) ; 
      break;

    default:
      WNSPARSESOLVE__ERROR( "No such solver" );

  }

  /*============================================================================
  // Allocate memory for parameters and workspace.
  //==========================================================================*/

  if( self->iDebug == 1 ) {
    printf( "Size of workspace = %d\n\n", a_ipar[3] );
  }

  /*============================================================================
  // Print matrix information, if desired.
  //==========================================================================*/

  if( self->iDebug == 1 ) {
    printf( "Number of rows in matrix = %lu\n", (unsigned long) i_rows );
    printf(
      "Number of non-zero elements in matrix = %lu\n\n",
      (unsigned long) i_elements
    );
  }

  /*============================================================================
  // Get work vectors.
  //==========================================================================*/

  if(
     (
       a_wk =
       ( double * ) malloc( ( (size_t) a_ipar[3] ) * sizeof( double ) )
     )
     == NULL
  )
  {
    WNSPARSESOLVE__ERROR( "Couldn't allocate memory in WnSparseSolve" );
  }

  p_in =
    gsl_vector_calloc( p_rhs_vector->size );

  /*============================================================================
  // Initialize arrays and iteration variables.
  //==========================================================================*/

  p_solution_vector =
    gsl_vector_alloc( p_rhs_vector->size );

  gsl_vector_memcpy( p_solution_vector, p_guess_vector );

  a_ipar[0] = 0;

  /*============================================================================
  // Loop.
  //==========================================================================*/

  while( a_ipar[6] <= self->iIterMax ) {

    switch( self->iSolverMethod )
    {

      case CG:
        cg_(
          &i_rows,
          p_rhs_vector->data,
          p_solution_vector->data,
          a_ipar,
          a_fpar,
          a_wk
        );
        break;

      case CGNR:
        cgnr_(
          &i_rows,
          p_rhs_vector->data,
          p_solution_vector->data,
          a_ipar,
          a_fpar,
          a_wk
        );
        break;

      case 2:
        bcg_(
          &i_rows,
          p_rhs_vector->data,
          p_solution_vector->data,
          a_ipar,
          a_fpar,
          a_wk
        );
        break;

      case DBCG:
        dbcg_(
          &i_rows,
          p_rhs_vector->data,
          p_solution_vector->data,
          a_ipar,
          a_fpar,
          a_wk
        );
        break;

      case BCGSTAB:
        bcgstab_(
          &i_rows,
          p_rhs_vector->data,
          p_solution_vector->data,
          a_ipar,
          a_fpar,
          a_wk
        );
        break;

      case TFQMR:
        tfqmr_(
          &i_rows,
          p_rhs_vector->data,
          p_solution_vector->data,
          a_ipar,
          a_fpar,
          a_wk
        );
        break;

      case FOM:
        fom_(
          &i_rows,
          p_rhs_vector->data,
          p_solution_vector->data,
          a_ipar,
          a_fpar,
          a_wk
        );
        break;

      case GMRES:
        gmres_(
          &i_rows,
          p_rhs_vector->data,
          p_solution_vector->data,
          a_ipar,
          a_fpar,
          a_wk
        );
        break;

      case FGMRES:
        fgmres_(
          &i_rows,
          p_rhs_vector->data,
          p_solution_vector->data,
          a_ipar,
          a_fpar,
          a_wk
        );
        break;

      case DQGMRES:
        dqgmres_(
          &i_rows,
          p_rhs_vector->data,
          p_solution_vector->data,
          a_ipar,
          a_fpar,
          a_wk
        );
        break;

      default:
        WNSPARSESOLVE__ERROR( "No solver found" );
        break;

    }

    if( i_iter == 0 ) {
      p_tmp =
        WnMatrix__computeMatrixTimesVector( p_matrix, p_solution_vector );
      gsl_vector_sub( p_tmp, p_rhs_vector );
      gsl_blas_ddot( p_tmp, p_tmp, &d_residual ); 
      gsl_vector_free( p_tmp );
      d_residual = sqrt( d_residual );
    }
    else {
      d_residual = a_fpar[4];
    } 

    if( i_iter != a_ipar[6] ) {
      if( self->iDebug == 1 && a_ipar[2] != 999 ) {
        fprintf(
          stdout, 
          "Iteration = %d   Residual = %e\n",
          a_ipar[6],
          d_residual
        );
        i_iter = a_ipar[6];
      }
    }

  /*---------------------------------------------------------------------------
  // Update in vector.
  //--------------------------------------------------------------------------*/

    if( a_ipar[7] > 0 ) {
      for( i = 0; i < p_in->size ; i++ )
        gsl_vector_set( p_in, i, a_wk[(size_t) a_ipar[7] - 1 + i] );
    }

  /*---------------------------------------------------------------------------
  // Matrix multiply, preconditioner solve,  or convergence check.
  //--------------------------------------------------------------------------*/

    switch( a_ipar[0] ) {

      case 1:
        p_out = 
          WnMatrix__computeMatrixTimesVector(
            p_matrix, p_in
          );
        for( i = 0; i < i_rows; i++ )
          a_wk[(size_t) a_ipar[8] - 1 + i] = gsl_vector_get( p_out, i );
        gsl_vector_free( p_out );
        break;

      case 2:
        p_out =
          WnMatrix__computeTransposeMatrixTimesVector(
            p_matrix, p_in
          );
        for( i = 0; i < i_rows; i++ )
          a_wk[(size_t) a_ipar[8] - 1 + i] = gsl_vector_get( p_out, i );
        gsl_vector_free( p_out );
        break;

      case 3: case 5:
        p_out =
          self->pfPreconditionerSolve(
            p_in,
            self->pPreconditionerUserData
        );
        for( i = 0; i < i_rows; i++ )
          a_wk[(size_t) a_ipar[8] - 1 + i] = gsl_vector_get( p_out, i );
        gsl_vector_free( p_out );
        break;

      case 4: case 6:
        p_out =
          self->pfPreconditionerTransposeSolve(
            p_in,
            self->pPreconditionerUserData
          );
        for( i = 0; i < i_rows; i++ )
          a_wk[(size_t) a_ipar[8] - 1 + i] = gsl_vector_get( p_out, i );
        gsl_vector_free( p_out );
        break;

      case 10:
        p_delta = gsl_vector_alloc( p_solution_vector->size );
        for( i = 0; i < p_solution_vector->size; i++ )
          gsl_vector_set( p_delta, i, a_wk[(size_t) a_ipar[8] - 1 + i ] );
        p_working_solution = gsl_vector_alloc( p_solution_vector->size );
        for( i = 0; i < p_solution_vector->size; i++ )
          gsl_vector_set(
            p_working_solution,
            i,
            a_wk[(size_t) a_ipar[7] - 1 + i]
          );

        a_ipar[10] = 
           self->pfConvergenceTester(
             p_delta,
             p_working_solution,
             self->pConvergenceTesterUserData
           );

        gsl_vector_free( p_delta );
        gsl_vector_free( p_working_solution );
        break;

      default:
        switch( a_ipar[0] )
        {

          case 0:
            if( self->iDebug == 1 ) {
              fprintf(
                stdout,
                "\nIterative solver has satisfied convergence test\n"
              );
              if( a_ipar[2] != 999 ) {
                fprintf(
                  stdout, 
                  "Target residual was %e\n",
                  a_fpar[3]
                );
              }
            }
            free( a_wk );   
            gsl_vector_free( p_in );
            return p_solution_vector;

          case -1:
            if( self->iDebug == 1 ) {
              fprintf(
                stdout,
                "\nIterative solver has iterated too many times.\n"
              );
              fprintf(
                stdout,
                "Maximum number of iterations = %d\n",
                self->iIterMax
              );
              fprintf(
                stdout,
                "Number of iterations performed = %d\n",
                a_ipar[6]
              );
            }
            free( a_wk );
            gsl_vector_free( p_in );
            gsl_vector_free( p_solution_vector );
            return NULL;

          case -2:
            if( self->iDebug == 1 ) {
              fprintf(
                stdout,
                "\nIterative solver was not given enough work space.\n"
              );
              fprintf(
                stdout,
                "The work space should have at least have %d elements.\n",
                a_ipar[3]
              );
            }
            free( a_wk );
            gsl_vector_free( p_in );
            gsl_vector_free( p_solution_vector );
            return NULL;

          case -3:
            if( self->iDebug == 1 ) {
              fprintf(
                stdout,
                "\nIterative solver is facing a breakdown.\n"
              );
            }
            free( a_wk );
            gsl_vector_free( p_in );
            gsl_vector_free( p_solution_vector );
            return NULL;

          case -4:
            if( self->iDebug == 1 ) {
              fprintf(
                stdout,
                "\nBoth tolerances are <= 0.\n"
              );
            }
            free( a_wk );
            gsl_vector_free( p_in );
            gsl_vector_free( p_solution_vector );
            return NULL;

          default:
            if( self->iDebug == 1 ) {
              fprintf(
                stdout,
                "\nIterative solve terminated with code = %d\n",
                a_ipar[0]
              );
            }
            free( a_wk );
            gsl_vector_free( p_in );
            gsl_vector_free( p_solution_vector );
            return NULL;

        }


    }

  }

  if( self->iDebug == 1 ) {
    fprintf(
      stdout,
      "\nIterative solver has iterated too many times.\n"
    );
    fprintf(
      stdout,
      "Maximum number of iterations = %d\n",
      self->iIterMax
    );
    fprintf(
      stdout,
      "Number of iterations performed = %d\n",
      a_ipar[6]
    );
  }
  free( a_wk );
  gsl_vector_free( p_in );
  gsl_vector_free( p_solution_vector );
  return NULL;

}

/*##############################################################################
// WnSparseSolve__Mat__free().
//############################################################################*/

void
WnSparseSolve__Mat__free(
  WnSparseSolve__Mat *self
)
{

  free( self );

}

/*##############################################################################
// WnSparseSolve__Mat__defaultSolve().
//############################################################################*/

gsl_vector *
WnSparseSolve__Mat__defaultSolve(
  gsl_vector *p_input_vector,
  void *p_user_data
)
{

  gsl_vector *p_output_vector;

  if( p_user_data )
    WNSPARSESOLVE__ERROR(
      "No user data should be passed into default solver.\n"
    );

  p_output_vector = gsl_vector_alloc( p_input_vector->size );

  gsl_vector_memcpy( p_output_vector, p_input_vector );

  return p_output_vector;

}

/*##############################################################################
// WnSparseSolve__Mat__defaultTransposeSolve().
//############################################################################*/

gsl_vector *
WnSparseSolve__Mat__defaultTransposeSolve(
  gsl_vector *p_input_vector,
  void *p_user_data
)
{

  gsl_vector *p_output_vector;

  if( p_user_data )
    WNSPARSESOLVE__ERROR(
      "No user data should be passed into default solver.\n"
    );

  p_output_vector = gsl_vector_alloc( p_input_vector->size );

  gsl_vector_memcpy( p_output_vector, p_input_vector );

  return p_output_vector;

}

/*##############################################################################
// WnSparseSolve__Exp__new().
//############################################################################*/

WnSparseSolve__Exp *
WnSparseSolve__Exp__new( void )
{

  WnSparseSolve__Exp *self;

  self = ( WnSparseSolve__Exp * ) malloc( sizeof( WnSparseSolve__Exp ) );

  self->iWorkSpace = EXP_WORKSPACE;
  self->iIterMax = EXP_ITER_MAX;
  self->iDebug = DEBUG_INIT;
  self->dTolerance = EXP_TOL;

  return self;

}
  
/*##############################################################################
// WnSparseSolve__Exp__setDebug().
//############################################################################*/

void
WnSparseSolve__Exp__setDebug(
  WnSparseSolve__Exp *self
)
{

  if( !self ) 
    WNSPARSESOLVE__ERROR( "Invalid input" );

  self->iDebug = 1;

}

/*##############################################################################
// WnSparseSolve__Exp__clearDebug().
//############################################################################*/

void
WnSparseSolve__Exp__clearDebug(
  WnSparseSolve__Exp *self
)
{

  if( !self ) 
    WNSPARSESOLVE__ERROR( "Invalid input" );

  self->iDebug = 0;

}

/*##############################################################################
// WnSparseSolve__Exp__updateMaximumIteration().
//############################################################################*/

void
WnSparseSolve__Exp__updateMaximumIterations(
  WnSparseSolve__Exp *self, int i_iter_max
)
{

  if( !self ) 
    WNSPARSESOLVE__ERROR( "Invalid input" );

  if( i_iter_max <= 0 )
    WNSPARSESOLVE__ERROR( "Maximum number of iterations must be > 0" );

  self->iIterMax = i_iter_max;

}

/*##############################################################################
// WnSparseSolve__Exp__updateWorkSpace().
//############################################################################*/

void
WnSparseSolve__Exp__updateWorkSpace(
  WnSparseSolve__Exp *self, int i_workspace
)
{

  if( !self ) 
    WNSPARSESOLVE__ERROR( "Invalid input" );

  if( i_workspace <= 0 )
    WNSPARSESOLVE__ERROR( "Workspace must be > 0" );

  if( i_workspace > 60 )
    WNSPARSESOLVE__ERROR( "Workspace must be <= 60" );

  self->iWorkSpace = i_workspace;

}

/*##############################################################################
// WnSparseSolve__Exp__updateTolerance().
//############################################################################*/

void
WnSparseSolve__Exp__updateTolerance(
  WnSparseSolve__Exp *self, double d_tolerance
)
{

  if( !self ) 
    WNSPARSESOLVE__ERROR( "Invalid input" );

  if( d_tolerance < 0. )
    WNSPARSESOLVE__ERROR( "Tolerance must be non-negative" );

  self->dTolerance = d_tolerance;

}

/*##############################################################################
// WnSparseSolve__Exp__free().
//############################################################################*/

void
WnSparseSolve__Exp__free(
  WnSparseSolve__Exp *self
)
{

  free( self );

}

/*##############################################################################
// WnSparseSolve__Exp__solve().
//############################################################################*/

gsl_vector *
WnSparseSolve__Exp__solve(
  WnSparseSolve__Exp *self,
  WnMatrix *p_matrix,
  const gsl_vector *p_input_vector,
  double d_tn
)
{

  size_t i_rows, n_fpar = 6;
  int i_indic, i_ierr, i_iter, i_workspace;
  int i_ipar[5] = {0,0,0,0,0};
  double d_tol;
  gsl_vector *p_x, *p_y, *p_u, *p_output_vector, *p_fpar;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if( !p_matrix || !p_input_vector )
    WNSPARSESOLVE__ERROR( "Invalid input" );

  /*============================================================================
  // Find debugging information.
  //==========================================================================*/

  i_ipar[3] = self->iDebug;

  if( i_ipar[3] == 1 ) printf( "\nUsing exppro\n\n" );

  /*============================================================================
  // Get number of rows.
  //==========================================================================*/

  i_rows = WnMatrix__getNumberOfRows( p_matrix );

  /*============================================================================
  // Set local variables from structure.
  //==========================================================================*/

  i_workspace = self->iWorkSpace;
  d_tol = self->dTolerance;

  /*============================================================================
  // Allocate memory for vectors.
  //==========================================================================*/

  p_x = gsl_vector_calloc( i_rows );
  p_y = gsl_vector_calloc( i_rows );
  p_u = gsl_vector_calloc( i_rows * ( (size_t) i_workspace + 1 ) );
  p_fpar = gsl_vector_calloc( n_fpar );
  p_output_vector = gsl_vector_alloc( p_input_vector->size );

  if( !p_x || !p_y || !p_u || !p_output_vector || !p_fpar )
      WNSPARSESOLVE__ERROR( "Couldn't allocate memory in WnSparseSolve" );

  gsl_vector_memcpy( p_output_vector, p_input_vector );

  i_indic = 0;
  i_iter = 0;

  while( i_indic != 1 ) {

    if( i_iter == self->iIterMax ) {
      if( self->iDebug == 1 )
        fprintf( stdout, "Exceeded maximum number of iterations.\n" );
      gsl_vector_free( p_x );
      gsl_vector_free( p_y );
      gsl_vector_free( p_u );
      gsl_vector_free( p_fpar );
      gsl_vector_free( p_output_vector );
      return NULL;
    }

    exppro_(
      &i_rows,
      &i_workspace,
      &d_tol,
      &d_tn,
      p_u->data,
      p_output_vector->data,
      p_x->data,
      p_y->data,
      p_fpar->data,
      i_ipar,
      &i_indic,
      &i_ierr
    );

    gsl_vector_free( p_y );

    p_y =
      WnMatrix__computeMatrixTimesVector(
        p_matrix,
        p_x
      );

    gsl_vector_scale( p_y, -1. );

    ++i_iter;

  }

  if( self->iDebug == 1 )
    printf( "%d iterations.\n", i_iter );

  /*============================================================================
  // Clean up and return.
  //==========================================================================*/

  gsl_vector_free( p_x );
  gsl_vector_free( p_y );
  gsl_vector_free( p_u );
  gsl_vector_free( p_fpar );

  return p_output_vector;

}

/*##############################################################################
// WnSparseSolve__Phi__new().
//############################################################################*/

WnSparseSolve__Phi *
WnSparseSolve__Phi__new( void )
{

  WnSparseSolve__Phi *self;

  self = ( WnSparseSolve__Phi * ) malloc( sizeof( WnSparseSolve__Phi ) );

  self->iWorkSpace = PHI_WORKSPACE;
  self->iIterMax = PHI_ITER_MAX;
  self->iDebug = DEBUG_INIT;
  self->dTolerance = PHI_TOL;

  return self;

}
  
/*##############################################################################
// WnSparseSolve__Phi__setDebug().
//############################################################################*/

void
WnSparseSolve__Phi__setDebug(
  WnSparseSolve__Phi *self
)
{

  if( !self ) 
    WNSPARSESOLVE__ERROR( "Invalid input" );

  self->iDebug = 1;

}

/*##############################################################################
// WnSparseSolve__Phi__clearDebug().
//############################################################################*/

void
WnSparseSolve__Phi__clearDebug(
  WnSparseSolve__Phi *self
)
{

  if( !self ) 
    WNSPARSESOLVE__ERROR( "Invalid input" );

  self->iDebug = 0;

}

/*##############################################################################
// WnSparseSolve__Phi__free().
//############################################################################*/

void
WnSparseSolve__Phi__free(
  WnSparseSolve__Phi *self
)
{

  free( self );

}

/*##############################################################################
// WnSparseSolve__Phi__updateMaximumIteration().
//############################################################################*/

void
WnSparseSolve__Phi__updateMaximumIterations(
  WnSparseSolve__Phi *self, int i_iter_max
)
{

  if( !self ) 
    WNSPARSESOLVE__ERROR( "Invalid input" );

  if( i_iter_max <= 0 )
    WNSPARSESOLVE__ERROR( "Maximum number of iterations must be > 0" );

  self->iIterMax = i_iter_max;

}

/*##############################################################################
// WnSparseSolve__Phi__updateWorkSpace().
//############################################################################*/

void
WnSparseSolve__Phi__updateWorkSpace(
  WnSparseSolve__Phi *self, int i_workspace
)
{

  if( !self ) 
    WNSPARSESOLVE__ERROR( "Invalid input" );

  if( i_workspace <= 0 )
    WNSPARSESOLVE__ERROR( "Workspace must be > 0" );

  if( i_workspace > 60 )
    WNSPARSESOLVE__ERROR( "Workspace must be <= 60" );

  self->iWorkSpace = i_workspace;

}

/*##############################################################################
// WnSparseSolve__Phi__updateTolerance().
//############################################################################*/

void
WnSparseSolve__Phi__updateTolerance(
  WnSparseSolve__Phi *self, double d_tolerance
)
{

  if( !self ) 
    WNSPARSESOLVE__ERROR( "Invalid input" );

  if( d_tolerance < 0. )
    WNSPARSESOLVE__ERROR( "Tolerance must be non-negative" );

  self->dTolerance = d_tolerance;

}

/*##############################################################################
// WnSparseSolve__Phi__solve().
//############################################################################*/

gsl_vector *
WnSparseSolve__Phi__solve(
  WnSparseSolve__Phi *self,
  WnMatrix *p_matrix,
  const gsl_vector *p_input_vector,
  const gsl_vector *p_const_vector,
  double d_tn
)
{

  size_t i_rows, n_fpar = 6;
  int i_indic, i_ierr, i_iter, i_workspace;
  int i_ipar[5] = {0,0,0,0,0};
  double d_tol;
  gsl_vector *p_x, *p_y, *p_u, *p_output_vector, *p_fpar;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if( !p_matrix || !p_input_vector || !p_const_vector )
    WNSPARSESOLVE__ERROR( "Invalid input" );

  /*============================================================================
  // Find debugging information.
  //==========================================================================*/

  if( self->iDebug == 1 )
      printf( "\nUsing phipro\n\n" );

  /*============================================================================
  // Set debugging flag.
  //==========================================================================*/

  i_ipar[3] = self->iDebug;

  /*============================================================================
  // Get number of rows.
  //==========================================================================*/

  i_rows = WnMatrix__getNumberOfRows( p_matrix );

  /*============================================================================
  // Set local variables from structure.
  //==========================================================================*/

  i_workspace = self->iWorkSpace;
  d_tol = self->dTolerance;

  /*============================================================================
  // Allocate memory for vectors.
  //==========================================================================*/

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

    p_y =
      WnMatrix__computeMatrixTimesVector(
        p_matrix,
        p_x
      );

    gsl_vector_scale( p_y, -1. );

    ++i_iter;

  }

  if( self->iDebug == 1 )
    printf( "%d iterations.\n", i_iter );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  gsl_vector_free( p_x );
  gsl_vector_free( p_y );
  gsl_vector_free( p_u );
  gsl_vector_free( p_fpar );

  return p_output_vector;

}

