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
//       Example to demonstrate how to use SPARSKIT incomplete lu
//       decomposition solvers with wn_sparse_solve.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include "ilu_solvers.h"

/*##############################################################################
// ilu_preconditioner_new().
//############################################################################*/

ilu_structure *
ilu_preconditioner_new( WnMatrix *p_matrix )
{

  size_t i, i_rows, i_w, i_lfil;
  int i_err, *a_jw;
  double *a_w;
  WnMatrix__Csr *p_csr;
  double d_drptol = D_DROP_TOL;
  int *a_row, *a_col;
  ilu_structure *self;

  /*============================================================================
  // Get CSR.
  //==========================================================================*/

  p_csr = WnMatrix__getCsr( p_matrix );

  /*============================================================================
  // Do allocations.
  //==========================================================================*/

  i_rows = WnMatrix__getNumberOfRows( p_matrix );

  i_lfil = WnMatrix__getNumberOfElements( p_matrix ) / i_rows + I_DELTA;

  i_w = 2 * ( i_rows + 1 + i_rows * i_lfil );

  self = ( ilu_structure * ) malloc( sizeof( ilu_structure ) );

  self->aA0 =
    ( double * )
    calloc(
      sizeof( double ), i_w
    );

  self->aJa0 = ( int * ) calloc( sizeof( int ), i_w );

  self->aIa0 =
    ( int * ) calloc(
      sizeof( int ), ( i_rows + 1 )
    );

  a_w =
    ( double * )
    calloc( sizeof( double ), i_rows + 1 );

  a_jw =
    ( int * )
    calloc( sizeof( int ), ( i_rows * 2 ) );

  /*============================================================================
  // Allocate row pointer and column arrays.  Do this for compatibility with
  // 64 bit compilers because of int and int mismatch.
  //==========================================================================*/

  a_row =
    ( int * )
    calloc( sizeof( long ), ( i_rows + 1 ) );

  a_col =
    ( int * )
    calloc( sizeof( long ), ( WnMatrix__getNumberOfElements( p_matrix ) ) );

  /*============================================================================
  // Assign row pointer and column arrays.
  //==========================================================================*/

  for( i = 0; i <= i_rows; i++ )
    a_row[i] = (int) WnMatrix__Csr__getRowPointerVector( p_csr )[i] + 1;

  for( i = 0; i < WnMatrix__getNumberOfElements( p_matrix ); i++ )
    a_col[i] = (int) WnMatrix__Csr__getColumnVector( p_csr )[i];

  /*============================================================================
  // Call pre-conditioner set up.
  //==========================================================================*/

  ilut_(
    &i_rows,
    WnMatrix__Csr__getValueVector( p_csr ),
    a_col,
    a_row,
    &i_lfil,
    &d_drptol,
    self->aA0,
    self->aJa0,
    self->aIa0,
    &i_w,
    a_w,
    a_jw,
    &i_err
  );

  free( a_row );
  free( a_col );
  free( a_w );
  free( a_jw );

  WnMatrix__Csr__free( p_csr );

  return self;

}

/*##############################################################################
// ilu_preconditioner_free().
//############################################################################*/

void
ilu_preconditioner_free( ilu_structure *self )
{

  free( self->aA0 );
  free( self->aJa0 );
  free( self->aIa0 );
  free( self );

}

/*##############################################################################
// my_ilu_solver().
//############################################################################*/

gsl_vector *
my_ilu_solver(
  gsl_vector *p_input_vector,
  ilu_structure *p_user_data
)
{

  size_t i_rows;
  gsl_vector *p_output_vector;

  i_rows = WnMatrix__get_gsl_vector_size( p_input_vector );

  p_output_vector = gsl_vector_alloc( i_rows );

  lusol_(
    &i_rows,
    WnMatrix__get_gsl_vector_array( p_input_vector ),
    WnMatrix__get_gsl_vector_array( p_output_vector ),
    p_user_data->aA0,
    p_user_data->aJa0,
    p_user_data->aIa0
  );

  return p_output_vector;

}

/*##############################################################################
// my_ilu_transpose_solver().
//############################################################################*/

gsl_vector *
my_ilu_transpose_solver(
  gsl_vector *p_input_vector,
  ilu_structure *p_user_data
)
{

  size_t i_rows;
  gsl_vector *p_output_vector;

  i_rows = WnMatrix__get_gsl_vector_size( p_input_vector );

  p_output_vector = gsl_vector_alloc( i_rows );

  lutsol_(
    &i_rows,
    WnMatrix__get_gsl_vector_array( p_input_vector ),
    WnMatrix__get_gsl_vector_array( p_output_vector ),
    p_user_data->aA0,
    p_user_data->aJa0,
    p_user_data->aIa0
  );

  return p_output_vector;

}
