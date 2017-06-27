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
//
//   <description>
//     <abstract>
//       Example of a user-supplied Coulomb correction factor function.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/
 
#include "coul_corr.h"

/*##############################################################################
// User-supplied Coulomb correction factor function.  Libnucnet uses GNU GSL
// definitions for physical constants, so the user-supplied routine should
// as well for consistency.
//############################################################################*/

double
my_coulomb_correction(
  Libnucnet__Species *p_species,
  double d_t9,
  double d_rho,
  double d_ye,
  user_coul_corr_data *p_data
) {

  double d_A3;
  double d_a_e, d_Gamma_e, d_Gamma_i;
  double d_term1, d_term2, d_term2a, d_term2b, d_term3a, d_term3b, d_term3;

  if( !p_species )
  {
    fprintf( stderr, "Invalid species!\n" );
    exit( EXIT_FAILURE );
  }

  if( !p_data )
  {
    fprintf( stderr, "Missing extra data for my_coulomb_correction().\n" );
    exit( EXIT_FAILURE );
  }

  d_A3   =  -(sqrt(3.0)/2.0) - ( p_data->dA1 / sqrt( p_data->dA2 ) );

  d_a_e = 4.0 * M_PI * d_rho * GSL_CONST_NUM_AVOGADRO * d_ye / 3.0; 
  d_a_e = pow( d_a_e, 1.0 / 3.0 );

  d_Gamma_e =
    gsl_pow_2(
      GSL_CONST_CGSM_ELECTRON_CHARGE * GSL_CONST_CGSM_SPEED_OF_LIGHT / 10.
    ) * d_a_e /
    ( GSL_CONST_CGSM_BOLTZMANN * d_t9 * GSL_CONST_NUM_GIGA );

  d_Gamma_i =
    d_Gamma_e *
    pow( ( double ) Libnucnet__Species__getZ( p_species ), 5.0 / 3.0);

  d_term1 =
    sqrt(
      d_Gamma_i*( p_data->dA2 + d_Gamma_i )
    );

  d_term2a = sqrt( d_Gamma_i / p_data->dA2 );
  d_term2b = sqrt( 1.0 + d_Gamma_i / p_data->dA2 );
  d_term2  = ( p_data->dA2 ) * log( d_term2a + d_term2b );

  d_term3a = sqrt( d_Gamma_i );
  d_term3b =
    atan(
      sqrt( d_Gamma_i )
    );
  d_term3 = 2.0*d_A3*( d_term3a - d_term3b );

  return ( p_data->dA1 ) *( d_term1 - d_term2 + d_term3 );

}

