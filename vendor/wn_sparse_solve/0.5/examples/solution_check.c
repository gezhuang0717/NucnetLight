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
//       Routine to print out and check solutions for matrix equations for
//       the wn_sparse_solve examples.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/
 
int
is_valid_solution( gsl_vector *, const char * );

int
is_valid_solution( gsl_vector *p_solution_vector, const char *s_debug )
{

  /*============================================================================
  // Check that the solution succeeded.
  //==========================================================================*/

  if( !p_solution_vector ) {
    if( !s_debug ) {
      fprintf(
        stderr,
        "\nSolution failed.  Try running with debug flag to find out why.\n\n"
      );
    }
    else if( strcmp( s_debug, "debug" ) == 0 ) {
      fprintf(
        stderr,
        "\nSolution failed.\n\n"
      );
    }
    else {
      fprintf(
        stderr,
        "\nSolution failed.  Try running with debug flag to find out why.\n\n"
      );
    }

    return 0;

  }

  return 1;
 
}
