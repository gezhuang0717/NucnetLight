/*//////////////////////////////////////////////////////////////////////////////
// <file type="public">
//   <license>
//     This file was originally written by George C. Jordan, IV and
//     Bradley S. Meyer.
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
//       Example of a user-supplied screening factor function
//       header file.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

/*##############################################################################
// Includes.
//############################################################################*/

#include <Libnucnet.h>

/*##############################################################################
// Use extern "C" for C++ compilers.
//############################################################################*/

#ifdef __cplusplus
extern "C"
{
#endif

/*##############################################################################
// Defines.
//############################################################################*/

#define D_THETA_E   0.0   /* setting to 1.0 for non-degenerate matter.
                             Set to 0.0 for the completely degenerate case.
                          */

/*##############################################################################
// User-supplied data structure.
//############################################################################*/

typedef struct user_screening_data {
  double dYe2;
} user_screening_data;

/*##############################################################################
// Prototypes for user-supplied functions.
//############################################################################*/

double
my_screening_function(
  double, double, double, int, int, int, int, user_screening_data *
);

double calculate_gamma_e( double, double, double );

double calculate_gamma_effective( int, int, double );

double
strong_screening_factor(
  int, int, int, int, double, double, double
);

double
weak_screening_factor(
  int, int, double, double, double, double, double
);

double intermediate_screening_factor( double, double );

#ifdef __cplusplus
} /* extern "C" */
#endif
