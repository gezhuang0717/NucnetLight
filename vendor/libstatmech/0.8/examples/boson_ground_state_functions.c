/*//////////////////////////////////////////////////////////////////////////////
// <file type="public">
//   <license>
//     This file was originally written by Tianhong Yu and Bradley S. Meyer.
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
//       File giving example ground-state boson functions.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libstatmech.h>
#include "boson_ground_state_functions.h"

/*##############################################################################
// boson_number_density_function()
//############################################################################*/

double
boson_number_density_function(
  Libstatmech__Boson *p_boson,
  double d_T,
  double d_alpha,
  double *p_volume
)
{

  if( d_T <= 0. )
  {
    fprintf( stderr, "Invalid input.\n" );
    exit( EXIT_FAILURE );
  }

  return
    Libstatmech__Boson__getMultiplicity( p_boson ) /
    (
      gsl_expm1( -d_alpha ) *
      *p_volume
    );

}

/*##############################################################################
// boson_pressure_function()
//############################################################################*/

double
boson_pressure_function(
  Libstatmech__Boson *p_boson,
  double d_T,
  double d_alpha,
  double *p_volume
)
{

  return
    -GSL_CONST_CGSM_BOLTZMANN * d_T *
    Libstatmech__Boson__getMultiplicity( p_boson) *
    log( -gsl_expm1( d_alpha ) ) /
    *p_volume;

}

/*##############################################################################
// boson_entropy_density_function() 
//############################################################################*/

double
boson_entropy_density_function(
  Libstatmech__Boson *p_boson,
  double d_T,
  double d_alpha,
  double *p_volume
)
{

  double d_result;

  if( !p_boson || d_T <= 0. )
  {
    fprintf( stderr, "Invalid input.\n" );
    exit( EXIT_FAILURE );
  }

  d_result =
    GSL_CONST_CGSM_BOLTZMANN *
    log( 
      -1 /
      gsl_expm1( d_alpha )
    ) 
    /
    *p_volume;
  
  d_result +=
    GSL_CONST_CGSM_BOLTZMANN *
    d_alpha /
    (  
      -gsl_expm1( -d_alpha ) *
      *p_volume
    );
  
  return d_result;

}

