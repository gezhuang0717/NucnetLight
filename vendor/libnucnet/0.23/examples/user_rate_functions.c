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
//       Example to demonstrate user-defined rate functions.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libnucnet__Reac.h>
#include "user_rate_functions.h"

/*##############################################################################
// register_my_rate_functions().
//############################################################################*/

void
register_my_rate_functions( Libnucnet__Reac *p_reac )
{

  /*============================================================================
  // Register Kunz et al. (2002) fit.
  //==========================================================================*/

  Libnucnet__Reac__registerUserRateFunction(
    p_reac,
    KUNZ_FIT,
    (Libnucnet__Reaction__userRateFunction) kunz_fit_function
  );

  /*============================================================================
  // Register Caughlan and Fowler (1988) weak fit.
  //==========================================================================*/

  Libnucnet__Reac__registerUserRateFunction(
    p_reac,
    CF88_WEAK_FIT,
    (Libnucnet__Reaction__userRateFunction) cf88_weak_fit_function
  );

  /*============================================================================
  // Register Caughlan and Fowler (1988) carbon fusion fit.
  //==========================================================================*/

  Libnucnet__Reac__registerUserRateFunction(
    p_reac,
    CF88_CARBON_FUSION_FIT,
    (Libnucnet__Reaction__userRateFunction) cf88_carbon_fusion_fit_function
  );

}

/*##############################################################################
// set_user_data_deallocators().
//############################################################################*/

int
set_user_data_deallocators(
  Libnucnet__Reac * p_reac
)
{

  return
    Libnucnet__Reac__setUserRateFunctionDataDeallocator(
      p_reac,
      CF88_WEAK_FIT,
      (Libnucnet__Reaction__user_rate_function_data_deallocator) free
    );

}

/*##############################################################################
// kunz_fit_function()
//############################################################################*/

double
kunz_fit_function(
  Libnucnet__Reaction *p_reaction,
  double d_t9,
  void *p_data
)
{

  double d_result;
  kunz_data k_data;

  if( p_data )
  {
    fprintf( stderr, "No extra data to this function.\n" );
    exit( EXIT_FAILURE );
  }

  k_data.iCount = 0;

  Libnucnet__Reaction__iterateUserRateFunctionProperties(
    p_reaction,
    "a",
    NULL,
    NULL,
    (Libnucnet__Reaction__user_rate_property_iterate_function)
      get_kunz_array,
    &k_data
  );

  if( k_data.iCount != 12 )
  {
    fprintf(
      stderr,
      "\nFor reaction %s:\n",
      Libnucnet__Reaction__getString( p_reaction )
    );
    fprintf(
      stderr,
      "  Invalid number of array elements in kunz_fit_function.\n\n"
    );
    exit( EXIT_FAILURE );
  }

  d_result =
    k_data.dA[0] /
    (
      gsl_pow_2( d_t9 ) *
      gsl_pow_2(
        1. + k_data.dA[1] * pow( d_t9, -2./3. )
      )
    )
    *
    exp(
      -k_data.dA[2] / pow( d_t9, 1./3. ) 
      - gsl_pow_2( d_t9 / k_data.dA[3]  )
    )
    +
    k_data.dA[4] * exp( -k_data.dA[6] / pow( d_t9, 1./3. ) )
    /
    (
      gsl_pow_2( d_t9 ) *
      gsl_pow_2(
        1. + k_data.dA[5] * pow( d_t9, -2./3. )
      )
    )
    +
    k_data.dA[7] * exp( -k_data.dA[8] / d_t9 ) / pow( d_t9, 3./2. )
    +
    k_data.dA[9] *
      (
        1. + k_data.dA[10] * pow( d_t9, 1./3. )
      ) /
      pow( d_t9, 2./3. )
    *
    exp( -k_data.dA[11] / pow( d_t9, 1./3. ) );

  return d_result;

}
    
/*##############################################################################
// get_kunz_array()
//############################################################################*/

void
get_kunz_array(
  const char *s_name,
  const char *s_tag1,
  const char *s_tag2,
  const char *s_property,
  kunz_data *p_kunz_data
)
{

  if( !s_name || s_tag2 )
  {
    fprintf( stderr, "Invalid input to get_kunz_array.\n" );
    exit( EXIT_FAILURE );
  }

  (*p_kunz_data).dA[atoi(s_tag1)] = atof( s_property );

  (*p_kunz_data).iCount++;

}

/*##############################################################################
// cf88_weak_fit_function()
//############################################################################*/

double
cf88_weak_fit_function(
  Libnucnet__Reaction *p_reaction,
  double d_t9,
  void *p_data
)
{

  double d_rate = 1., d_part = 0.;
  const char *s_value;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if( !p_reaction )
  {
    fprintf( stderr, "No reaction supplied.\n" );
    exit( EXIT_FAILURE );
  }

  if( d_t9 <= 0. )
  {
    fprintf( stderr, "Invalid temperature.\n" );
    exit( EXIT_FAILURE );
  }

  if( !p_data )
  {
    fprintf( stderr, "cf88_weak_fit_function needs extra data.\n" );
    exit( EXIT_FAILURE );
  }

  /*============================================================================
  // T9^(1./3.).
  //==========================================================================*/

  if( 
      (
        s_value =
          Libnucnet__Reaction__getUserRateFunctionProperty(
            p_reaction,
            "t913",
            NULL,
            NULL
          )
      )
  )
    d_rate += atof( s_value ) * pow( d_t9, 1./3. );

  /*============================================================================
  // T9^(2./3.).
  //==========================================================================*/

  if( 
      (
        s_value =
          Libnucnet__Reaction__getUserRateFunctionProperty(
            p_reaction,
            "t923",
            NULL,
            NULL
          )
      )
  )
    d_rate += atof( s_value ) * pow( d_t9, 2./3. );

  /*============================================================================
  // Exp.
  //==========================================================================*/

  if( 
      (
        s_value =
          Libnucnet__Reaction__getUserRateFunctionProperty(
            p_reaction,
            "t9m",
            NULL,
            NULL
          )
      )
  )
    d_part = atof( s_value ) / d_t9;

  if( 
      (
        s_value =
          Libnucnet__Reaction__getUserRateFunctionProperty(
            p_reaction,
            "t9expm",
            NULL,
            NULL
          )
      )
  )
    d_part *= exp( atof( s_value ) / d_t9 );

  d_rate += d_part;

  /*============================================================================
  // Prefactor.
  //==========================================================================*/

  if( 
      (
        s_value =
          Libnucnet__Reaction__getUserRateFunctionProperty(
            p_reaction,
            "t912m",
            NULL,
            NULL
          )
      )
  )
    d_rate *= atof( s_value ) / sqrt( d_t9 );

  /*============================================================================
  // Check rate.
  //==========================================================================*/

  d_rate *= *(double *) p_data;

  s_value =
    Libnucnet__Reaction__getUserRateFunctionProperty(
      p_reaction,
      "rate_max",
      NULL,
      NULL
    );

  if( s_value )
    d_rate = GSL_MIN_DBL( d_rate, atof( s_value ) );

  /*============================================================================
  // Done.
  //==========================================================================*/

  return d_rate;

}

/*##############################################################################
// cf88_carbon_fusion_fit_function().
//############################################################################*/

double
cf88_carbon_fusion_fit_function(
  Libnucnet__Reaction *p_reaction,
  double d_t9,
  void *p_data
)
{

  double d_rate = 0.;
  const char *s_fac = NULL;

  if( p_data )
  {
    fprintf( stderr, "No extra data for carbon fusion function.\n" );
    exit( EXIT_FAILURE );
  }

  d_rate = cf88_carbon_fusion( d_t9 );

  if( d_t9 >= 6. )
    s_fac =
      Libnucnet__Reaction__getUserRateFunctionProperty( 
        p_reaction,
        S1,
        NULL,
        NULL
      );
  else if( d_t9 < 6. && d_t9 > 0.0 )
    s_fac =
      Libnucnet__Reaction__getUserRateFunctionProperty( 
        p_reaction,
        S2,
        NULL,
        NULL
      );

  if( !s_fac )
  {
    fprintf( stderr, "Invalid factor in cf88_carbon_fusion_fit_function.\n" );
    exit( EXIT_FAILURE );
  }

  d_rate *= atof( s_fac );

  return d_rate;

}

/*##############################################################################
// cf88_carbon_fusion().
//############################################################################*/

double
cf88_carbon_fusion( double d_t9 )
{

  double d_t9a;

  d_t9a = d_t9 / ( 1. + 0.0369 * d_t9 );

  return
    4.27e26 *
    pow( d_t9a, 5./6. ) *
    exp( -84.165 / pow( d_t9a, 1./3. ) - 2.12e-3 * gsl_pow_3( d_t9 ) ) /
    pow( d_t9, 3./2. );

}
