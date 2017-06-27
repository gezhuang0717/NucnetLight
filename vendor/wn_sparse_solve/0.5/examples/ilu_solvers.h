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
//       Header file for use with ilu_solvers.c.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <WnSparseSolve.h>

#ifdef __cplusplus
extern "C" {
#endif

/*##############################################################################
// Defines.  User can change the values as appropriate.
//############################################################################*/

#define    D_DROP_TOL    0.     /* The drop tolerance. */
#define    I_DELTA       1      /* Extra amount to add to lfil. */ 

/*##############################################################################
// Data structure.  This passes the ilut preconditioner arrays into the solver
// routines.
//############################################################################*/

typedef struct ilu_structure {
  int *aIa0;
  int *aJa0;
  double *aA0;
} ilu_structure;

/*##############################################################################
// Prototypes.
//############################################################################*/

/*==============================================================================
// ilut() is the SPARSKIT2 routine that sets the incomplete lu decomposition.
//============================================================================*/

void
ilut_(
  size_t *,
  double *,
  int *,
  int *,
  size_t *,
  double *,
  double *,
  int *,
  int *,
  size_t *,
  double *,
  int *,
  int *
);

/*==============================================================================
// lusol() is the SPARSKIT2 routine that solves the matrix equation
// with the ilu preconditioner.
//============================================================================*/

void
lusol_(
  size_t *,
  double *,
  double *,
  double *,
  int *,
  int *
);

/*==============================================================================
// lusolt() is the SPARSKIT2 routine that solves the transpose matrix equation
// with the ilu preconditioner.
//============================================================================*/

void
lutsol_(
  size_t *,
  double *,
  double *,
  double *,
  int *,
  int *
);

/*==============================================================================
// ilu_preconditioner_new() creates a new ilu_structure from the matrix
// data.
//============================================================================*/

ilu_structure *
ilu_preconditioner_new( WnMatrix * );

/*==============================================================================
// ilu_preconditioner_free() deallocates the memory for an
// ilu_structure.
//============================================================================*/

void
ilu_preconditioner_free( ilu_structure * );

/*==============================================================================
// my_ilu_solver() calls lusol with the preconditioner data in
// ilu_structure.
//============================================================================*/

gsl_vector *
my_ilu_solver(
  gsl_vector *,
  ilu_structure *
);

/*==============================================================================
// my_ilu_transpose_solver() calls lutsol with the preconditioner data
// in ilu_structure.
//============================================================================*/

gsl_vector *
my_ilu_transpose_solver(
  gsl_vector *,
  ilu_structure *
);

#ifdef __cplusplus
}
#endif

