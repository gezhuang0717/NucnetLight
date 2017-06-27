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
//        Header file for the wn_sparse_solve source code.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#ifndef WN_SPARSE_SOLVE_H
#define WN_SPARSE_SOLVE_H

/*##############################################################################
// Use extern "C" for C++ compilers.
//############################################################################*/

#ifdef __cplusplus
extern "C"
{
#endif

/*##############################################################################
// Includes.
//############################################################################*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "WnMatrix.h"

#include <gsl/gsl_blas.h>

/*##############################################################################
// WnSparseSolve defines.
//############################################################################*/

/* Uncomment following line to debug */
/* #define DEBUG */

#ifndef DEBUG

#define WNSPARSESOLVE__ERROR( s_reason ) \
       do { \
         fprintf( stderr, "ERROR: %s in %s on line %d\n", \
           s_reason, __FILE__, __LINE__ \
         ); \
         exit( EXIT_FAILURE ); \
       } while (0)

#else

#define WNSPARSESOLVE__ERROR( s_reason ) \
       do { \
         fprintf( stderr, "ERROR: %s in %s on line %d\n", \
           s_reason, __FILE__, __LINE__ \
         ); \
         abort(); \
       } while (0)

#endif

enum {
  CG, CGNR, BCG, DBCG, BCGSTAB, TFQMR, FOM, GMRES, FGMRES, DQGMRES
};

enum {
  CONVERGE_INITIAL, CONVERGE_RHS
};

#define WN_SPARSE_SOLVE_BUF_SIZE  256

#define S_CG       "cg"
#define S_CGNR     "cgnr"
#define S_BCG      "bcg"
#define S_DBCG     "dbcg"
#define S_BCGSTAB  "bcgstab"
#define S_TFQMR    "tfqmr"
#define S_FOM      "fom"
#define S_GMRES    "gmres"
#define S_FGMRES   "fgmres"
#define S_DQGMRES  "dqgmres"

#define S_CONVERGE_INITIAL  "initial residual"
#define S_CONVERGE_RHS      "rhs"

#define DEFAULT_SOLVER   BCG
#define DEFAULT_CONVERGE CONVERGE_INITIAL
#define MAT_ITER_MAX     100
#define D_ABS_TOL        1.e-8
#define D_REL_TOL        1.e-8
#define DEBUG_INIT       0

#define EXP_TOL          1.e-8
#define EXP_ITER_MAX     100
#define EXP_WORKSPACE    40

#define PHI_TOL          1.e-8
#define PHI_ITER_MAX     100
#define PHI_WORKSPACE    40

/*##############################################################################
// Sparsekit prototypes.
//############################################################################*/

void cg_( size_t *, double [], double [], int [], double [], double * );
void cgnr_( size_t *, double [], double [], int [], double [], double * );
void bcg_( size_t *, double [], double [], int [], double [], double * );
void dbcg_( size_t *, double [], double [], int [], double [], double * );
void bcgstab_(
  size_t *, double [], double [], int [], double [], double *
);
void tfqmr_(
  size_t *, double [], double [], int [], double [], double *
);
void fom_( size_t *, double [], double [], int [], double [], double * );
void gmres_(
  size_t *, double [], double [], int [], double [], double *
);
void fgmres_(
  size_t *, double [], double [], int [], double [], double *
);
void dqgmres_(
  size_t *, double [], double [], int [], double [], double *
);

/*##############################################################################
// <user_routine name="WnSparseSolve__Mat__PreconditionerSolver()">
//
//   <description>
//     <abstract>
//       Optional user-supplied routine to compute the solution to the matrix
//       equation Mx = b, where the inverse of M is an approximation to
//       the inverse of the matrix for which the solution is required.
//     </abstract>
//     <keywords>
//       preconditioner, solver, user, supplied
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// gsl_vector *
// WnSparseSolve__Mat__PreconditionerSolver(
//   gsl_vector *p_rhs_vector,
//   void *p_user_data
// );
//     </calling_sequence>
//
//     <param
//       name="p_rhs_vector"
//       kind="in,positional,required"
//     >
//       A pointer to the gsl_vector containing the right-hand-side vector
//       in the matrix equation Mx = b.
//     </param>
//
//     <param
//       name="p_user_data"
//       kind="in,positional,required"
//     >
//       A pointer to the user's data structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       User's routine must return a new gsl_vector giving the solution to
//       Mx = b, where the inverse of M is the approximate
//       inverse matrix (M or its inverse is passed in through the user's data
//       structure) and b is the right-hand-side vector.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//   </usage>
//
// </user_routine>
//############################################################################*/

typedef gsl_vector *( *WnSparseSolve__Mat__PreconditionerSolver ) (
  gsl_vector *, void *
);

/*##############################################################################
// <user_routine name="WnSparseSolve__Mat__PreconditionerTransposeSolver()">
//
//   <description>
//     <abstract>
//       Optional user-supplied routine to compute the solution to the matrix
//       equation (transpose of M)x = b, where the inverse of M is an
//       approximation to the inverse of the matrix for which the solution is
//       required.
//     </abstract>
//     <keywords>
//       preconditioner, solver, user, supplied
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// gsl_vector *
// WnSparseSolve__Mat__PreconditionerTransposeSolver(
//   gsl_vector *p_rhs_vector,
//   void *p_user_data
// );
//     </calling_sequence>
//
//     <param
//       name="p_rhs_vector"
//       kind="in,positional,required"
//     >
//       A pointer to the gsl_vector containing the right-hand-side vector
//       in the matrix equation (transpose of M)x = b.
//     </param>
//
//     <param
//       name="p_user_data"
//       kind="in,positional,required"
//     >
//       A pointer to the user's data structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       User's routine must return a new gsl_vector giving the solution to
//       (transpose of M)x = b, where the inverse of M is the approximate
//       inverse matrix (M or its inverse is passed in through the user's
//       data structure) and b is the right-hand-side vector.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//   </usage>
//
// </user_routine>
//############################################################################*/

typedef gsl_vector *( *WnSparseSolve__Mat__PreconditionerTransposeSolver ) (
  gsl_vector *, void *
);

/*##############################################################################
// <user_routine name="WnSparseSolve__Mat__ConvergenceTester()">
//
//   <description>
//     <abstract>
//       Optional user-supplied routine to test the convergence of a solution.
//     </abstract>
//     <keywords>
//       preconditioner, solver, user, supplied
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// int
// WnSparseSolve__Mat__ConvergenceTester(
//   gsl_vector *p_change_vector,
//   gsl_vector *p_solution_vector,
//   void *p_user_data
// );
//     </calling_sequence>
//
//     <param
//       name="p_change_vector"
//       kind="in,positional,required"
//     >
//       A pointer to the gsl_vector containing the update to the solution
//       in an iteration.
//     </param>
//
//     <param
//       name="p_solution_vector"
//       kind="in,positional,required"
//     >
//       A pointer to the gsl_vector containing the current solution
//       in an iteration.
//     </param>
//
//     <param
//       name="p_user_data"
//       kind="in,positional,required"
//     >
//       A pointer to the user's data structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       User's routine must return a zero (0) if the solution does
//       not satisfy the user-defined convergence criterion or one (1)
//       if it does.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//   </usage>
//
// </user_routine>
//############################################################################*/

typedef int ( *WnSparseSolve__Mat__ConvergenceTester ) (
  gsl_vector *, gsl_vector *, void *
);

/*##############################################################################
// <class name="WnSparseSolve__Mat">
//
//   <description>
//     <abstract>
//       WnSparseSolve__Mat is a structure for handling settings for SPARSKIT2
//       iterative solvers.
//     </abstract>
//     <keywords>
//       sparse, matrix, webnucleo, solver
//     </keywords>
//   </description>
//
//   <authors>
//     <current>
//       <author userid="mbradle" start_date="2008/06/05"/>
//     </current>
//   </authors>
//
//   <compatibility>
//     Compiled successfully in gcc (GCC) 3.4.0.
//   </compatibility>
//
// </class>
//############################################################################*/

typedef struct WnSparseSolve__Mat {
  int iSolverMethod;
  int iConvergenceMethod;
  int iDebug;
  int iIterMax;
  double dRelTolerance;
  double dAbsTolerance;
  WnSparseSolve__Mat__PreconditionerSolver pfPreconditionerSolve;
  WnSparseSolve__Mat__PreconditionerTransposeSolver
    pfPreconditionerTransposeSolve;
  WnSparseSolve__Mat__ConvergenceTester
    pfConvergenceTester;
  void *pPreconditionerUserData;
  void *pConvergenceTesterUserData;
} WnSparseSolve__Mat;

/*##############################################################################
// <class name="WnSparseSolve__Exp">
//
//   <description>
//     <abstract>
//       WnSparseSolve__Exp is a structure for handling settings for SPARSKIT2
//       matrix exponentiation.
//     </abstract>
//     <keywords>
//       sparse, matrix, webnucleo, solver, exponentiation
//     </keywords>
//   </description>
//
//   <authors>
//     <current>
//       <author userid="mbradle" start_date="2008/06/05"/>
//     </current>
//   </authors>
//
//   <compatibility>
//     Compiled successfully in gcc (GCC) 3.4.0.
//   </compatibility>
//
// </class>
//############################################################################*/

typedef struct WnSparseSolve__Exp {
  int iDebug;
  int iIterMax;
  double dTolerance;
  int iWorkSpace;
} WnSparseSolve__Exp;

/*##############################################################################
// <class name="WnSparseSolve__Phi">
//
//   <description>
//     <abstract>
//       WnSparseSolve__Phi is a structure for handling settings for SPARSKIT2
//       matrix exponentiation with an additional constant added vector.
//     </abstract>
//     <keywords>
//       sparse, matrix, webnucleo, solver, exponentiation
//     </keywords>
//   </description>
//
//   <authors>
//     <current>
//       <author userid="mbradle" start_date="2008/06/05"/>
//     </current>
//   </authors>
//
//   <compatibility>
//     Compiled successfully in gcc (GCC) 3.4.0.
//   </compatibility>
//
// </class>
//############################################################################*/

typedef struct WnSparseSolve__Phi {
  int iDebug;
  int iIterMax;
  double dTolerance;
  int iWorkSpace;
} WnSparseSolve__Phi;

/*##############################################################################
// API routines.
//############################################################################*/

/*##############################################################################
// <routine name="WnSparseSolve__Mat__new()">
//
//   <description>
//     <abstract>
//       Create a new sparse solver.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, solver
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// WnSparseSolve__Mat *
// WnSparseSolve__Mat__new( );
//     </calling_sequence>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Upon successful return, a new solver has been created.  The default
//       values for the solver are that the type is a biconjugate gradient
//       solver, the maximum number of iterations allowed is 100,
//       the absolute and relative tolerances are 1.e-8, the
//       preconditioner matrix is simply the identity matrix, and debugging
//       is turned off.  These values can be updated by other API routines.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Create a new solver:
//       </synopsis>
//       <code>
// p_my_solver = WnSparseSolve__Mat__new();
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

WnSparseSolve__Mat * WnSparseSolve__Mat__new( void );

/*##############################################################################
// <routine name="WnSparseSolve__Mat__solve()">
//
//   <description>
//     <abstract>
//       Solve a matrix equation.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// gsl_vector *
// WnSparseSolve__Mat__solve(
//   WnSparseSolve__Mat *self,
//   WnMatrix *p_matrix,
//   gsl_vector *p_rhs,
//   gsl_vector *p_guess
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="self"
//     >
//       A pointer to a WnSparseSolve__Mat structure.
//     </param>
//     <param
//       name="p_matrix"
//       kind="in,positional,required"
//       doc="matrix"
//     >
//       A pointer to a WnMatrix structure containing the matrix data for
//       the matrix equation.
//     </param>
//     <param
//       name="p_rhs"
//       kind="in,positional,required"
//       doc="in_vector"
//     >
//       A gsl_vector giving the right hand side vector in the matrix
//       equation.
//     </param>
//     <param
//       name="p_guess"
//       kind="in,positional,required"
//       doc="in_vector"
//     >
//       A gsl_vector giving a guess vector for the solution to the matrix
//       equation.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="pre" id="in_vector">
//        The vector should have a number of elements equal to the number
//        of columns of the matrix.
//     </doc>
//
//     <doc kind="post" id="result">
//       Routine returns a gsl_vector with the solution (if solution successful)
//       or NULL otherwise.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Solve the matrix equation A x = b where A is stored as
//          WnMatrix *p_my_matrix, b is stored in the gsl_vector p_rhs, and
//          a guess for the solution is stored in p_guess.  Use the
//          solver WnSparseSolve__Mat *p_my_solver and store
//          result in the gsl_vector p_sol:
//       </synopsis>
//       <code>
// p_sol =
//   WnSparseSolve__Mat__solve(
//     p_my_solver, p_my_matrix, p_rhs, p_guess
//   );
//
// if( p_sol )
//     fprintf( stdout, "Solution succeeded!\n" );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

gsl_vector *
WnSparseSolve__Mat__solve(
  WnSparseSolve__Mat *,
  WnMatrix *,
  gsl_vector *,
  gsl_vector *
);

/*##############################################################################
// <routine name="WnSparseSolve__Mat__setDebug()">
//
//   <description>
//     <abstract>
//       Set the debug flag for a sparse matrix solver.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, debug
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnSparseSolve__Mat__setDebug(
//   WnSparseSolve__Mat *self
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="solver"
//     >
//       A pointer to a WnSparseSolve__Mat structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Upon successful return, the debug flag has been set for the solver
//       so that debugging information will be printed to stdout during
//       solution iterations.  If the input is invalid, error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Set the debug flag for WnSparseSolve__Mat *p_my_solver:
//       </synopsis>
//       <code>
//   WnSparseSolve__Mat__setDebug(
//     p_my_solver
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnSparseSolve__Mat__setDebug( WnSparseSolve__Mat * );

/*##############################################################################
// <routine name="WnSparseSolve__Mat__clearDebug()">
//
//   <description>
//     <abstract>
//       Clear the debug flag for a sparse matrix solver.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, debug
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnSparseSolve__Mat__clearDebug(
//   WnSparseSolve__Mat *self
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="solver"
//     >
//       A pointer to a WnSparseSolve__Mat structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Upon successful return, the debug flag has been cleared for the solver
//       so that debugging information will no longer be printed to stdout
//       during solution iterations.  If the input is invalid, error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Clear the debug flag for WnSparseSolve__Mat *p_my_solver:
//       </synopsis>
//       <code>
//   WnSparseSolve__Mat__clearDebug(
//     p_my_solver
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnSparseSolve__Mat__clearDebug( WnSparseSolve__Mat * );

/*##############################################################################
// <routine name="WnSparseSolve__Mat__updateSolverMethod()">
//
//   <description>
//     <abstract>
//       Update the solver method.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, method
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnSparseSolve__Mat__updateSolverMethod(
//   WnSparseSolve__Mat *self,
//   const char *s_method
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="solver"
//     >
//       A pointer to a WnSparseSolve__Mat structure.
//     </param>
//
//     <param
//       name="s_method"
//       kind="in,positional,required"
//       doc="method"
//     >
//       A string giving the name of the solver method to use.  The possible
//       values are 1) "cg"--conjugate gradient, 2) "cgnr"--conjugate gradient
//       with normalized residue, 3) "bcg"--bi-conjugate gradient,
//       4) "dbcg"--bi-conjugate gradient with partial pivoting,
//       5) "bcgstab"--stabilized bi-conjugate gradient,
//       6) "tfqmr"--transpose-free quasi-minimum residual, 7) "fom"--full
//       orthogonalization, 8) "gmres"--generalized minimum residual,
//       9) "fgmres"--flexible generalized minimum residual, or
//       10) "dqgmres"--direct version of quasi generalized minimum residual.
//       See the SPARSKIT2 documentation for more details.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Upon successful return, the solver method has been updated
//       to the input value.  If the input is invalid, error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Set the solver method for WnSparseSolve__Mat *p_my_solver to
//          the full orthogonalization method:
//       </synopsis>
//       <code>
//   WnSparseSolve__Mat__updateSolverMethod(
//     p_my_solver, "fom"
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnSparseSolve__Mat__updateSolverMethod(
  WnSparseSolve__Mat *, const char *
);

/*##############################################################################
// <routine name="WnSparseSolve__Mat__updateConvergenceMethod()">
//
//   <description>
//     <abstract>
//       Update the convergence method.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, convergence, method
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnSparseSolve__Mat__updateConvergenceMethod(
//   WnSparseSolve__Mat *self,
//   const char *s_method
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="solver"
//     >
//       A pointer to a WnSparseSolve__Mat structure.
//     </param>
//
//     <param
//       name="s_method"
//       kind="in,positional,required"
//       doc="method"
//     >
//       A string giving the name of the convergence method to use.
//       The current possible values are 1) "initial residual"--the
//       solver iterates until ||residual|| <= rtol * ||initial residual||
//       + atol or 2) "rhs"--the solver iterates until ||residual||
//       <= rtol * ||rhs|| + atol.  Here ||.|| denotes the 2-norm, rtol
//       is the relative tolerance, and atol is the absolute tolerance.
//       This flag is ignored if a WnSparseSolve__Mat__ConvergenceTester
//       is specified with WnSparseSolve__Mat__updateConvergenceTester().
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Upon successful return, the solver method has been updated
//       to the input value.  If the input is invalid, error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Set the solver method for WnSparseSolve__Mat *p_my_solver to
//          use the rhs criterion:
//       </synopsis>
//       <code>
//   WnSparseSolve__Mat__updateConvergenceMethod(
//     p_my_solver, "rhs"
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnSparseSolve__Mat__updateConvergenceMethod(
  WnSparseSolve__Mat *, const char *
);

/*##############################################################################
// <routine name="WnSparseSolve__Mat__getSolverMethod()">
//
//   <description>
//     <abstract>
//       Retrieve the current solver method.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, method
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// const char *
// WnSparseSolve__Mat__getSolverMethod(
//   WnSparseSolve__Mat *self
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="solver"
//     >
//       A pointer to a WnSparseSolve__Mat structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       The routine returns a string indicating the current solver method.
//       If the input is invalid, error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Print the current solver method for WnSparseSolve__Mat *p_my_solver:
//       </synopsis>
//       <code>
// printf(
//   "The current solver is %s\n",  
//   WnSparseSolve__Mat__getSolverMethod(
//     p_my_solver
//   )
// );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

const char *
WnSparseSolve__Mat__getSolverMethod(
  WnSparseSolve__Mat *
);

/*##############################################################################
// <routine name="WnSparseSolve__Mat__updateMaximumIterations()">
//
//   <description>
//     <abstract>
//       Update the maximum number of iterations allowed for a solver.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, relative, maximum, iterations
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnSparseSolve__Mat__updateMaximumIterations(
//   WnSparseSolve__Mat *self,
//   int i_iter_max
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="solver"
//     >
//       A pointer to a WnSparseSolve__Mat structure.
//     </param>
//
//     <param
//       name="i_iter_max"
//       kind="in,positional,required"
//       doc="iter_max"
//     >
//       An int giving the new maximum number of iterations.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="pre" id="iter_max">
//       The maximum number of iterations must be greater than 0.
//     </doc>
//
//     <doc kind="post" id="result">
//       Upon successful return, the maximum number of iterations has been
//       updated to the input value. If the input is invalid, error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Set the maximum number of iterations for WnSparseSolve__Mat
//          *p_my_solver to 50:
//       </synopsis>
//       <code>
//   WnSparseSolve__Mat__updateMaximumIterations(
//     p_my_solver, 50
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnSparseSolve__Mat__updateMaximumIterations(
  WnSparseSolve__Mat *, int
);

/*##############################################################################
// <routine name="WnSparseSolve__Mat__updateRelativeTolerance()">
//
//   <description>
//     <abstract>
//       Update the relative tolerance for a solver.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, relative, tolerance
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnSparseSolve__Mat__updateRelativeTolerance
//   WnSparseSolve__Mat *self,
//   double d_relative_tolerance
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="solver"
//     >
//       A pointer to a WnSparseSolve__Mat structure.
//     </param>
//
//     <param
//       name="d_relative_tolerance"
//       kind="in,positional,required"
//       doc="rel_tolerance"
//     >
//       A double giving the new tolerance.
//     </param>
//
//     <doc kind="pre" id="rel_tolerance">
//       The tolerance should be a non-negative number.
//     </doc>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Upon successful return, the relative tolerance has been updated
//       to the input value. If the input is invalid, error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Set the relative tolerance for WnSparseSolve__Mat *p_my_solver
//          to 1.e-12:
//       </synopsis>
//       <code>
//   WnSparseSolve__Mat__updateRelativeTolerance(
//     p_my_solver, 1.e-12
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnSparseSolve__Mat__updateRelativeTolerance(
  WnSparseSolve__Mat *, double
);

/*##############################################################################
// <routine name="WnSparseSolve__Mat__updateAbsoluteTolerance()">
//
//   <description>
//     <abstract>
//       Update the relative tolerance for a solver.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, relative, tolerance
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnSparseSolve__Mat__updateAbsoluteTolerance
//   WnSparseSolve__Mat *self,
//   double d_relative_tolerance
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="solver"
//     >
//       A pointer to a WnSparseSolve__Mat structure.
//     </param>
//
//     <param
//       name="d_relative_tolerance"
//       kind="in,positional,required"
//       doc="abs_tolerance"
//     >
//       A double giving the new tolerance.
//     </param>
//
//     <doc kind="pre" id="abs_tolerance">
//       The tolerance should be a non-negative number.
//     </doc>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Upon successful return, the absolute tolerance has been updated
//       to the input value. If the input is invalid, error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Set the absolute tolerance for WnSparseSolve__Mat *p_my_solver
//          to 1.e-6:
//       </synopsis>
//       <code>
//   WnSparseSolve__Mat__updateAbsoluteTolerance(
//     p_my_solver, 1.e-6
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnSparseSolve__Mat__updateAbsoluteTolerance(
  WnSparseSolve__Mat *, double
);

/*##############################################################################
// <routine name="WnSparseSolve__Mat__updatePreconditionerSolver()">
//
//   <description>
//     <abstract>
//       Set the preconditioner solver.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, preconditioner
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnSparseSolve__Mat__updatePreconditionerSolver(
//   WnSparseSolve__Mat *self,
//   WnSparseSolve__Mat__PreconditionerSolver pf_preconditioner_solver
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="solver"
//     >
//       A pointer to a WnSparseSolve__Mat structure.
//     </param>
//
//     <param
//       name="pf_preconditioner_solver"
//       kind="in,positional,required"
//       doc="pre_solver"
//     >
//       The name of the user-supplied preconditioner solver.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Upon successful return, the preconditioner solver has been updated
//       to the input value. If the input is invalid, error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Set the preconditioner solver for WnSparseSolve__Mat *p_my_solver to
//          my_preconditioner_solver:
//       </synopsis>
//       <code>
//   WnSparseSolve__Mat__updatePreconditionerSolver(
//     p_my_solver,
//     (WnSparseSolve__Mat__PreconditionerSolver)
//        my_preconditioner_solver
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnSparseSolve__Mat__updatePreconditionerSolver(
  WnSparseSolve__Mat *,
  WnSparseSolve__Mat__PreconditionerSolver
);

/*##############################################################################
// <routine name="WnSparseSolve__Mat__updatePreconditionerTransposeSolver()">
//
//   <description>
//     <abstract>
//       Set the preconditioner transpose solver.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, preconditioner, transpose
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnSparseSolve__Mat__updatePreconditionerTransposeSolver(
//   WnSparseSolve__Mat *self,
//   WnSparseSolve__Mat__PreconditionerTransposeSolver
//     pf_preconditioner_transpose_solver
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="solver"
//     >
//       A pointer to a WnSparseSolve__Mat structure.
//     </param>
//
//     <param
//       name="pf_preconditioner_transpose_solver"
//       kind="in,positional,required"
//       doc="pre_solver"
//     >
//       The name of the user-supplied preconditioner transpose solver.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Upon successful return, the preconditioner transpose solver has
//       been updated to the input value. If the input is invalid,
//       error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Set the preconditioner transpose solver for
//          WnSparseSolve *p_my_solver to
//          my_preconditioner_transpose solver:
//       </synopsis>
//       <code>
//   WnSparseSolve__Mat__updatePreconditionerTransposeSolver(
//     p_my_solver,
//     (WnSparseSolve__Mat__PreconditionerTransposeSolver)
//        my_preconditioner_transpose_solver
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnSparseSolve__Mat__updatePreconditionerTransposeSolver(
  WnSparseSolve__Mat *,
  WnSparseSolve__Mat__PreconditionerTransposeSolver
);

/*##############################################################################
// <routine name="WnSparseSolve__Mat__updateConvergenceTester()">
//
//   <description>
//     <abstract>
//       Set the convergence tester.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, convergence, tester
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnSparseSolve__Mat__updateConvergenceTester(
//   WnSparseSolve__Mat *self,
//   WnSparseSolve__Mat__ConvergenceTester pf_convergence_tester
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="solver"
//     >
//       A pointer to a WnSparseSolve__Mat structure.
//     </param>
//
//     <param
//       name="pf_convergence_tester"
//       kind="in,positional,required"
//       doc="convergence"
//     >
//       The name of the user-supplied convergence tester.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Upon successful return, the convergence tester has been updated
//       to the input value. If the input is invalid, error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Set the convergence tester for WnSparseSolve__Mat *p_my_solver to
//          my_tester:
//       </synopsis>
//       <code>
//   WnSparseSolve__Mat__updateConvergenceTester(
//     p_my_solver,
//     (WnSparseSolve__Mat__ConvergenceTester)
//        my_tester
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnSparseSolve__Mat__updateConvergenceTester(
  WnSparseSolve__Mat *,
  WnSparseSolve__Mat__ConvergenceTester
);

/*##############################################################################
// <routine name="WnSparseSolve__Mat__updatePreconditionerUserData()">
//
//   <description>
//     <abstract>
//       Set the user-supplied data to accompany the user's preconditioner
//       solvers.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, preconditioner, transpose, user, data
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnSparseSolve__Mat__updatePreconditionerUserData(
//   WnSparseSolve__Mat *self,
//   void *p_user_data
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="solver"
//     >
//       A pointer to a WnSparseSolve__Mat structure.
//     </param>
//
//     <param
//       name="p_user_data"
//       kind="in,positional,required"
//       doc="pre_solver"
//     >
//       A pointer to the user-supplied data structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Upon successful return, the user data has
//       been updated to the input value. If the input is invalid,
//       error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Set the extra data to be supplied by the user (which is
//          stored in the data structure my_data) for
//          WnSparseSolve__Mat *p_my_solver:
//       </synopsis>
//       <code>
//   WnSparseSolve__Mat__updatePreconditionerUserData(
//     p_my_solver,
//     &my_data
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnSparseSolve__Mat__updatePreconditionerUserData(
  WnSparseSolve__Mat *, void *
);

/*##############################################################################
// <routine name="WnSparseSolve__Mat__updateConvergenceTesterUserData()">
//
//   <description>
//     <abstract>
//       Set the user-supplied data to accompany the user's convergence
//       tester.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, convergence, tester, user, data
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnSparseSolve__Mat__updateConvergenceTesterUserData(
//   WnSparseSolve__Mat *self,
//   void *p_user_data
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="solver"
//     >
//       A pointer to a WnSparseSolve__Mat structure.
//     </param>
//
//     <param
//       name="p_user_data"
//       kind="in,positional,required"
//       doc="pre_solver"
//     >
//       A pointer to the user-supplied data structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Upon successful return, the user data has
//       been updated to the input value. If the input is invalid,
//       error handling is invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Set the extra data to be supplied by the user for
//          the convergence tester (which is
//          stored in the data structure my_data) for
//          WnSparseSolve__Mat *p_my_solver:
//       </synopsis>
//       <code>
//   WnSparseSolve__Mat__updateConvergenceTesterUserData(
//     p_my_solver,
//     &my_data
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnSparseSolve__Mat__updateConvergenceTesterUserData(
  WnSparseSolve__Mat *, void *
);

/*##############################################################################
// <routine name="WnSparseSolve__Mat__free()">
//
//   <description>
//     <abstract>
//       Free the memory allocated for a WnSparseSolve__Mat structure.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, free
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnSparseSolve__Mat__free(
//   WnSparseSolve__Mat *self
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="solver"
//     >
//       A pointer to a WnSparseSolve__Mat structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Upon successful return, the solver has been freed. The solver
//       pointer self is no longer available for use--it must be reallocated
//       using WnSparseSolve__Mat__new().
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Free the memory for WnSparseSolve__Mat *p_my_solver:
//       </synopsis>
//       <code>
//   WnSparseSolve__Mat__free(
//     p_my_solver
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnSparseSolve__Mat__free(
  WnSparseSolve__Mat *
);

/*##############################################################################
// <routine name="WnSparseSolve__Exp__new()">
//
//   <description>
//     <abstract>
//       Create a new sparse exponentiation solver.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, solver, exponentiation
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// WnSparseSolve__Exp *
// WnSparseSolve__Exp__new( );
//     </calling_sequence>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Upon successful return, a new solver has been created.  The default
//       values for the solver are that the maximum number of iterations
//       allowed is 100, the tolerances are 1.e-8, the workspace has
//       size 40, and debugging is turned off.  These values can be updated
//       by other API routines.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Create a new exponentiation solver:
//       </synopsis>
//       <code>
// p_my_exp_solver = WnSparseSolve__Exp__new();
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

WnSparseSolve__Exp *
WnSparseSolve__Exp__new( void );

/*##############################################################################
// <routine name="WnSparseSolve__Exp__setDebug()">
//
//   <description>
//     <abstract>
//       Set the debug flag for a sparse matrix exponentiation solver.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, debug, exponentiation
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnSparseSolve__Exp__setDebug(
//   WnSparseSolve__Exp *self
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="solver"
//     >
//       A pointer to a WnSparseSolve__Exp structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Upon successful return, the debug flag has been set for the solver
//       so that debugging information will be printed to stdout during
//       solution iterations.  If the input is invalid, error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Set the debug flag for WnSparseSolve__Exp *p_my_exp_solver:
//       </synopsis>
//       <code>
//   WnSparseSolve__Exp__setDebug(
//     p_my_exp_solver
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnSparseSolve__Exp__setDebug( WnSparseSolve__Exp * );

/*##############################################################################
// <routine name="WnSparseSolve__Exp__clearDebug()">
//
//   <description>
//     <abstract>
//       Clear the debug flag for a sparse matrix exponentiation solver.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, debug, exponentiation
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnSparseSolve__Exp__clearDebug(
//   WnSparseSolve__Exp *self
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="solver"
//     >
//       A pointer to a WnSparseSolve__Exp structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Upon successful return, the debug flag has been cleared for the solver
//       so that debugging information will no longer be printed to stdout
//       during solution iterations.  If the input is invalid, error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Clear the debug flag for WnSparseSolve__Exp *p_my_exp_solver:
//       </synopsis>
//       <code>
//   WnSparseSolve__Exp__clearDebug(
//     p_my_exp_solver
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnSparseSolve__Exp__clearDebug( WnSparseSolve__Exp * );

/*##############################################################################
// <routine name="WnSparseSolve__Exp__updateMaximumIterations()">
//
//   <description>
//     <abstract>
//       Update the maximum number of iterations allowed for an exponentation
//        solver.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, relative, maximum, iterations,
//       exponentiation
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnSparseSolve__Exp__updateMaximumIterations(
//   WnSparseSolve__Exp *self,
//   int i_iter_max
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="solver"
//     >
//       A pointer to a WnSparseSolve__Exp structure.
//     </param>
//
//     <param
//       name="i_iter_max"
//       kind="in,positional,required"
//       doc="iter_max"
//     >
//       An int giving the new maximum number of iterations.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="pre" id="iter_max">
//       The maximum number of iterations must be greater than 0.
//     </doc>
//
//     <doc kind="post" id="result">
//       Upon successful return, the maximum number of iterations has been
//       updated to the input value. If the input is invalid, error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Set the maximum number of iterations for WnSparseSolve__Exp
//          *p_my_exp_solver to 50:
//       </synopsis>
//       <code>
//   WnSparseSolve__Mat__updateMaximumIterations(
//     p_my_exp_solver, 50
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnSparseSolve__Exp__updateMaximumIterations(
  WnSparseSolve__Exp *, int
);

/*##############################################################################
// <routine name="WnSparseSolve__Exp__updateWorkSpace()">
//
//   <description>
//     <abstract>
//       Update the allowed workspace for an exponentation solver.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, relative, maximum, workspace,
//       exponentiation
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnSparseSolve__Exp__updateWorkSpace(
//   WnSparseSolve__Exp *self,
//   int i_workspace
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="solver"
//     >
//       A pointer to a WnSparseSolve__Exp structure.
//     </param>
//
//     <param
//       name="i_workspace"
//       kind="in,positional,required"
//       doc="workspace"
//     >
//       An int giving the new workspace.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="pre" id="iter_max">
//       The workspace must be > 0 and <= 60.
//     </doc>
//
//     <doc kind="post" id="result">
//       Upon successful return, the workspace has been
//       updated to the input value. If the input is invalid, error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Set the workspace for WnSparseSolve__Exp *p_my_exp_solver to 10:
//       </synopsis>
//       <code>
//   WnSparseSolve__Mat__updateWorkSpace(
//     p_my_exp_solver, 10
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnSparseSolve__Exp__updateWorkSpace(
  WnSparseSolve__Exp *, int
);

/*##############################################################################
// <routine name="WnSparseSolve__Exp__updateTolerance()">
//
//   <description>
//     <abstract>
//       Update the tolerance for an exponentiation solver.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, exponentiation, tolerance
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnSparseSolve__Exp__updateTolerance
//   WnSparseSolve__Exp *self,
//   double d_tolerance
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="solver"
//     >
//       A pointer to a WnSparseSolve__Mat structure.
//     </param>
//
//     <param
//       name="d_tolerance"
//       kind="in,positional,required"
//       doc="tolerance"
//     >
//       A double giving the new tolerance.
//     </param>
//
//     <doc kind="pre" id="tolerance">
//       The tolerance should be a non-negative number.
//     </doc>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Upon successful return, the tolerance has been updated
//       to the input value. If the input is invalid, error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Set the tolerance for WnSparseSolve__Exp *p_my_exp_solver
//          to 1.e-12:
//       </synopsis>
//       <code>
//   WnSparseSolve__Exp__updateTolerance(
//     p_my_exp_solver, 1.e-12
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnSparseSolve__Exp__updateTolerance(
  WnSparseSolve__Exp *, double
);

/*##############################################################################
// <routine name="WnSparseSolve__Exp__free()">
//
//   <description>
//     <abstract>
//       Free the memory allocated for a WnSparseSolve__Exp structure.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, free
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnSparseSolve__Exp__free(
//   WnSparseSolve__Exp *self
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="solver"
//     >
//       A pointer to a WnSparseSolve__Exp structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Upon successful return, the solver has been freed. The solver
//       pointer self is no longer available for use--it must be reallocated
//       using WnSparseSolve__Exp__new().
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Free the memory for WnSparseSolve__Exp *p_my_exp_solver:
//       </synopsis>
//       <code>
//   WnSparseSolve__Exp__free(
//     p_my_exp_solver
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnSparseSolve__Exp__free( WnSparseSolve__Exp * );

/*##############################################################################
// <routine name="WnSparseSolve__Exp__solve()">
//
//   <description>
//     <abstract>
//       Solve a matrix equation dY/dt = AY by matrix exponentiation.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, exponentiation
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// gsl_vector *
// WnSparseSolve__Exp__solve(
//   WnSparseSolve__Exp *self,
//   WnMatrix *p_matrix,
//   gsl_vector *p_initial_vector,
//   double d_time
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="self"
//     >
//       A pointer to a WnSparseSolve__Exp structure.
//     </param>
//     <param
//       name="p_matrix"
//       kind="in,positional,required"
//       doc="matrix"
//     >
//       A pointer to a WnMatrix structure containing the matrix data 
//       (that is, the matrix A in the matrix equation dY/dt = AY).
//     </param>
//     <param
//       name="p_initial_vector"
//       kind="in,positional,required"
//       doc="in_vector"
//     >
//       A gsl_vector giving the initial vector in the matrix
//       equation.
//     </param>
//     <param
//       name="d_time"
//       kind="in,positional,required"
//       doc="d_time"
//     >
//       A double giving the time over which to integrate the equation.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="pre" id="in_vector">
//        The vector should have a number of elements equal to the number
//        of columns of the matrix.
//     </doc>
//
//     <doc kind="post" id="result">
//       Routine returns a gsl_vector with the solution (if solution successful)
//       or NULL otherwise.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Solve the matrix equation dY/dt = AY, where A is stored as
//          WnMatrix *p_my_matrix, over time 10.  The solution is
//          Y(t) = exp( A*10 ) Y(0).  Y(0) is the initial vector stored
//          in p_in.  Use the exponentiation solver p_my_exp and store
//          the result in the gsl_vector p_sol:
//       </synopsis>
//       <code>
// p_sol =
//   WnSparseSolve__Exp__solve(
//     p_my_exp, p_my_matrix, p_in, 10.
//   );
//
// if( p_sol )
//   fprintf( stdout, "Solution succeeded!\n" );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

gsl_vector *
WnSparseSolve__Exp__solve(
  WnSparseSolve__Exp *,
  WnMatrix *,
  const gsl_vector *,
  double
);

/*##############################################################################
// <routine name="WnSparseSolve__Phi__new()">
//
//   <description>
//     <abstract>
//       Create a new sparse exponentiation plus constant solver.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, solver, exponentiation
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// WnSparseSolve__Phi *
// WnSparseSolve__Phi__new( );
//     </calling_sequence>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Upon successful return, a new solver has been created.  The default
//       values for the solver are that the maximum number of iterations
//       allowed is 100, the tolerances are 1.e-8, the workspace has
//       size 40, and debugging is turned off.  These values can be updated
//       by other API routines.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//         Create a new exponentiation solver:
//       </synopsis>
//       <code>
// p_my_phi_solver = WnSparseSolve__Phi__new();
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

WnSparseSolve__Phi *
WnSparseSolve__Phi__new( void );

/*##############################################################################
// <routine name="WnSparseSolve__Phi__setDebug()">
//
//   <description>
//     <abstract>
//       Set the debug flag for a sparse matrix exponentiation plus
//       constant solver.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, debug, exponentiation
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnSparseSolve__Phi__setDebug(
//   WnSparseSolve__Phi *self
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="solver"
//     >
//       A pointer to a WnSparseSolve__Phi structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Upon successful return, the debug flag has been set for the solver
//       so that debugging information will be printed to stdout during
//       solution iterations.  If the input is invalid, error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Set the debug flag for WnSparseSolve__Phi *p_my_phi_solver:
//       </synopsis>
//       <code>
//   WnSparseSolve__Phi__setDebug(
//     p_my_phi_solver
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnSparseSolve__Phi__setDebug( WnSparseSolve__Phi * );

/*##############################################################################
// <routine name="WnSparseSolve__Phi__clearDebug()">
//
//   <description>
//     <abstract>
//       Clear the debug flag for a sparse matrix exponentiation plus
//       constant solver.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, debug, exponentiation
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnSparseSolve__Phi__clearDebug(
//   WnSparseSolve__Phi *self
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="solver"
//     >
//       A pointer to a WnSparseSolve__Phi structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Upon successful return, the debug flag has been cleared for the solver
//       so that debugging information will no longer be printed to stdout
//       during solution iterations.  If the input is invalid, error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Clear the debug flag for WnSparseSolve__Phi *p_my_phi_solver:
//       </synopsis>
//       <code>
//   WnSparseSolve__Phi__clearDebug(
//     p_my_phi_solver
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnSparseSolve__Phi__clearDebug( WnSparseSolve__Phi * );

/*##############################################################################
// <routine name="WnSparseSolve__Phi__updateMaximumIterations()">
//
//   <description>
//     <abstract>
//       Update the maximum number of iterations allowed for an exponentation
//        plus constant solver.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, relative, maximum, iterations,
//       exponentiation
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnSparseSolve__Phi__updateMaximumIterations(
//   WnSparseSolve__Phi *self,
//   int i_iter_max
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="solver"
//     >
//       A pointer to a WnSparseSolve__Phi structure.
//     </param>
//
//     <param
//       name="i_iter_max"
//       kind="in,positional,required"
//       doc="iter_max"
//     >
//       An int giving the new maximum number of iterations.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="pre" id="iter_max">
//       The maximum number of iterations must be greater than 0.
//     </doc>
//
//     <doc kind="post" id="result">
//       Upon successful return, the maximum number of iterations has been
//       updated to the input value. If the input is invalid, error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Set the maximum number of iterations for WnSparseSolve__Phi
//          *p_my_phi_solver to 50:
//       </synopsis>
//       <code>
//   WnSparseSolve__Mat__updateMaximumIterations(
//     p_my_phi_solver, 50
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnSparseSolve__Phi__updateMaximumIterations(
  WnSparseSolve__Phi *, int
);

/*##############################################################################
// <routine name="WnSparseSolve__Phi__updateWorkSpace()">
//
//   <description>
//     <abstract>
//       Update the allowed workspace for an exponentation plus constant solver.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, relative, maximum, workspace,
//       exponentiation
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnSparseSolve__Phi__updateWorkSpace(
//   WnSparseSolve__Phi *self,
//   int i_workspace
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="solver"
//     >
//       A pointer to a WnSparseSolve__Phi structure.
//     </param>
//
//     <param
//       name="i_workspace"
//       kind="in,positional,required"
//       doc="workspace"
//     >
//       An int giving the new workspace.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="pre" id="iter_max">
//       The workspace must be > 0 and <= 60.
//     </doc>
//
//     <doc kind="post" id="result">
//       Upon successful return, the workspace has been
//       updated to the input value. If the input is invalid, error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Set the workspace for WnSparseSolve__Phi *p_my_phi_solver to 10:
//       </synopsis>
//       <code>
//   WnSparseSolve__Mat__updateWorkSpace(
//     p_my_phi_solver, 10
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnSparseSolve__Phi__updateWorkSpace(
  WnSparseSolve__Phi *, int
);

/*##############################################################################
// <routine name="WnSparseSolve__Phi__updateTolerance()">
//
//   <description>
//     <abstract>
//       Update the tolerance for an exponentiation plus constant solver.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, exponentiation, tolerance
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnSparseSolve__Phi__updateTolerance
//   WnSparseSolve__Phi *self,
//   double d_tolerance
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="solver"
//     >
//       A pointer to a WnSparseSolve__Mat structure.
//     </param>
//
//     <param
//       name="d_tolerance"
//       kind="in,positional,required"
//       doc="tolerance"
//     >
//       A double giving the new tolerance.
//     </param>
//
//     <doc kind="pre" id="tolerance">
//       The tolerance should be a non-negative number.
//     </doc>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Upon successful return, the tolerance has been updated
//       to the input value. If the input is invalid, error handling is
//       invoked.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Set the tolerance for WnSparseSolve__Phi *p_my_phi_solver
//          to 1.e-12:
//       </synopsis>
//       <code>
//   WnSparseSolve__Phi__updateTolerance(
//     p_my_phi_solver, 1.e-12
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnSparseSolve__Phi__updateTolerance(
  WnSparseSolve__Phi *, double
);

/*##############################################################################
// <routine name="WnSparseSolve__Phi__free()">
//
//   <description>
//     <abstract>
//       Free the memory allocated for a WnSparseSolve__Phi structure.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, free
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// void
// WnSparseSolve__Phi__free(
//   WnSparseSolve__Phi *self
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="solver"
//     >
//       A pointer to a WnSparseSolve__Phi structure.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="post" id="result">
//       Upon successful return, the solver has been freed. The solver
//       pointer self is no longer available for use--it must be reallocated
//       using WnSparseSolve__Phi__new().
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Free the memory for WnSparseSolve__Phi *p_my_phi_solver:
//       </synopsis>
//       <code>
//   WnSparseSolve__Phi__free(
//     p_my_solver
//   );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

void
WnSparseSolve__Phi__free( WnSparseSolve__Phi * );

/*##############################################################################
// <routine name="WnSparseSolve__Phi__solve()">
//
//   <description>
//     <abstract>
//       Solve a matrix equation dY/dt = AY + P, where P is a constant vector,
//       by matrix exponentiation.
//     </abstract>
//     <keywords>
//       sparse, webnucleo, matrix, solve, exponentiation, constant vector
//     </keywords>
//   </description>
//
//   <usage>
//
//     <calling_sequence>
// gsl_vector *
// WnSparseSolve__Phi__solve(
//   WnSparseSolve__Phi *self,
//   WnMatrix *p_matrix,
//   gsl_vector *p_initial_vector,
//   gsl_vector *p_constant_vector
// );
//       
//     </calling_sequence>
//
//     <param
//       name="self"
//       kind="in,positional,required"
//       doc="self"
//     >
//       A pointer to a WnSparseSolve__Phi structure.
//     </param>
//     <param
//       name="p_matrix"
//       kind="in,positional,required"
//       doc="matrix"
//     >
//       A pointer to a WnMatrix structure containing the matrix data 
//       (that is, the matrix A in the equation dY/dt = AY + P).
//     </param>
//     <param
//       name="p_initial_vector"
//       kind="in,positional,required"
//       doc="in_vector"
//     >
//       A gsl_vector giving the initial vector in the matrix
//       equation.
//     </param>
//     <param
//       name="p_constant_vector"
//       kind="in,positional,required"
//       doc="const_vector"
//     >
//       A gsl_vector giving the constant vector in the matrix
//       equation.
//     </param>
//
//     <param
//       kind="return"
//       doc="result"
//     />
//
//     <doc kind="pre" id="in_vector">
//        The vector should have a number of elements equal to the number
//        of columns of the matrix.
//     </doc>
//
//     <doc kind="post" id="result">
//       Routine returns a gsl_vector with the solution (if solution successful)
//       or NULL otherwise.
//     </doc>
//
//     <doc kind="bug" >
//     </doc>
//
//     <doc kind="example" id="example">
//       <synopsis>
//          Solve the matrix equation dY/dt = AY + P where A is stored as
//          WnMatrix *p_my_matrix over time 1.  Y(0) is the initial vector
//          stored in p_in.  The constant vector is p_constant.
//          Use the exponentiation solver p_my_phi and store
//          the result in the gsl_vector p_sol:
//       </synopsis>
//       <code>
// p_sol =
//   WnSparseSolve__Phi__solve(
//     p_my_phi, p_my_matrix, p_in, p_const, 1.
//   );
//
// if( p_sol )
//   fprintf( stdout, "Solution succeeded!\n" );
//       </code>
//     </doc>
//
//   </usage>
//
// </routine>
//############################################################################*/

gsl_vector *
WnSparseSolve__Phi__solve(
  WnSparseSolve__Phi *,
  WnMatrix *,
  const gsl_vector *,
  const gsl_vector *,
  double
);

/*#############################################################################/
// Non-API routines.
//############################################################################*/

gsl_vector *
WnSparseSolve__Mat__defaultSolve(
  gsl_vector *,
  void *
);

gsl_vector *
WnSparseSolve__Mat__defaultTransposeSolve(
  gsl_vector *,
  void *
);

/*#############################################################################/
// Declare Sparskit prototypes
//############################################################################*/

void exppro_(
  size_t *,
  int *,
  double *,
  double *,
  double[],
  double[],
  double[],
  double[],
  double[],
  int[],
  int *,
  int *
);
  
void phipro_(
  size_t *,
  int *,
  double *,
  double *,
  double[],
  double[],
  double[],
  double[],
  double[],
  double[],
  int[],
  int *,
  int *
);

#ifdef __cplusplus
}
#endif

#endif  /* WN_SPARSE_SOLVE_H */
