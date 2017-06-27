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
//       Example to demonstrate a user-defined convergence test.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include "my_convergence_tester.h"

/*##############################################################################
// my_convergence_tester().
//############################################################################*/

int
my_convergence_tester(
  gsl_vector *p_change,
  gsl_vector *p_solution,
  void *p_user_data
)
{

  size_t i;
  double d_checkT, d_check = 0;
  struct my_convergence_data *p_data =
    ( struct my_convergence_data * ) p_user_data;

  if( p_data->dCutOff < 0. ) {
    fprintf( stderr, "CutOff must >= 0!\n" );
    exit( EXIT_FAILURE );
  }

  for( i = 0; i < WnMatrix__get_gsl_vector_size( p_solution ); i++ )
  {
    if( fabs( gsl_vector_get( p_solution, i ) ) > p_data->dCutOff ) {
      d_checkT =
        fabs(
          gsl_vector_get( p_change, i ) /
          gsl_vector_get( p_solution, i )
        );
      if( d_checkT > d_check ) d_check = d_checkT;
    }
  }

  if( p_data->iDebug ) {
     fprintf(
       stdout,
       "Check = %e\n",
       d_check
     );
  }

  if( d_check < p_data->dConvergence )
    return 1;
  else
    return 0;

} 
