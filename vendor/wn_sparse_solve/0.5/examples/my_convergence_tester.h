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
//       Header file for use with my_convergence_tester.c.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <WnSparseSolve.h>

#ifdef __cplusplus
extern "C" {
#endif

/*##############################################################################
// User data structure.
//############################################################################*/

struct my_convergence_data{
  int iDebug;
  double dCutOff;
  double dConvergence;
};

/*##############################################################################
// prototype.
//############################################################################*/

int
my_convergence_tester(
  gsl_vector *, gsl_vector *, void *
);

#ifdef __cplusplus
}
#endif

