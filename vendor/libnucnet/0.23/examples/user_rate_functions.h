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
//       Header file for the user-defined rate functions.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libnucnet__Reac.h>

/*##############################################################################
// Use extern "C" for C++ compilers.
//############################################################################*/

#ifdef __cplusplus
extern "C"
{
#endif

/*##############################################################################
// Globals.
//############################################################################*/

typedef struct {
  double dA[12];
  int iCount;
} kunz_data;

/*##############################################################################
// Defines
//############################################################################*/

#define KUNZ_FIT                 "kunz fit"
#define CF88_WEAK_FIT            "cf88 weak fit"
#define CF88_CARBON_FUSION_FIT   "cf88 carbon fusion fit"
#define S1                       "f_t9_ge_6"
#define S2                       "f_3.3_le_t9_lt_6"

/*##############################################################################
// Prototypes.
//############################################################################*/

void
register_my_rate_functions( Libnucnet__Reac * );

int
set_user_data_deallocators( Libnucnet__Reac * );

double
kunz_fit_function(
  Libnucnet__Reaction *,
  double,
  void *
);

void
get_kunz_array(
  const char *,
  const char *,
  const char *,
  const char *,
  kunz_data *
);

double
cf88_weak_fit_function(
  Libnucnet__Reaction *,
  double,
  void *
);

double
cf88_carbon_fusion_fit_function(
  Libnucnet__Reaction *,
  double,
  void *
);

double
cf88_carbon_fusion( double );

#ifdef __cplusplus
} /* extern "C" */
#endif
