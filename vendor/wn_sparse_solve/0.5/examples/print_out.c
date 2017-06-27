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
//       Routine to print out matrix and vector information for
//       examples.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

/*##############################################################################
// Prototype.
//############################################################################*/

#ifdef __cplusplus
extern "C" {
#endif

void print_out( WnMatrix *, gsl_vector * );
void solution_print( WnMatrix *, gsl_vector *, gsl_vector * );

#ifdef __cplusplus
}
#endif

/*##############################################################################
// print_out().
//############################################################################*/

void print_out( WnMatrix *p_matrix, gsl_vector *p_vector )
{

  WnMatrix__Line *p_row;
  size_t i, j;

  /*============================================================================
  // Print out matrix.
  //==========================================================================*/

  fprintf( stdout, "\nMatrix:\n\n" );
  fprintf( stdout, "Row   Column       Value\n" );
  fprintf( stdout, "---   ------       -----\n" );

  for( i = 1; i <= WnMatrix__getNumberOfRows( p_matrix); i++ )
  {

    p_row = WnMatrix__getRow( p_matrix, i );

    for( j = 0; j < WnMatrix__Line__getNumberOfElements( p_row ); j++ )
    {
      fprintf(
        stdout,
        "%3lu%9lu%12.2e\n",
        (unsigned long) i,
        (unsigned long) WnMatrix__Line__getNonZeroIndices( p_row )[j],
        WnMatrix__Line__getNonZeroElements( p_row )[j]
      );
    }

    WnMatrix__Line__free( p_row );

  } 
      

  /*============================================================================
  // Print out vector.
  //==========================================================================*/

  fprintf( stdout, "\nRight-hand-side vector:\n\n" );
  for(
    i = 0;
    i < WnMatrix__get_gsl_vector_size( p_vector );
    i++
  ) {
    fprintf(
      stdout,
      "i = %d  rhs[i] = %g\n",
      (unsigned int) i,
      gsl_vector_get( p_vector, i )
    );
  }

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
      "i = %d  a_sol[i] = %15.10e\n",
      (unsigned int) i_rows,
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
      "i = %d  matrix * a_sol[i] = %f  rhs[i] = %f\n",
      (unsigned int) i_rows,
      gsl_vector_get( p_check_vector, i_rows ),
      gsl_vector_get( p_rhs_vector, i_rows )
    );
  }

  fprintf( stdout, "\n" );

  gsl_vector_free( p_check_vector );

}
