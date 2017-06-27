/*//////////////////////////////////////////////////////////////////////////////
// <file type="public">
//   <license>
//      Copyright (c) 2008-2012 Clemson University.
//
//      This file is part of the Clemson Webnucleo group's
//      libstatmech module, originally developed by Tianhong Yu
//      and Bradley S. Meyer.  For more information,
//      please see http://www.webnucleo.org.
//
//      This is free software; you can redistribute it and/or modify it
//      under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 2 of the License, or
//      (at your option) any later version.
//
//      This software is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//
//      You should have received a copy of the GNU General Public License
//      along with this software (please see the "gnu_gpl.txt" file in the doc/
//      directory of this distribution); if not, write to the Free Software
//      Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
//      USA
//   </license>
// </file>
//
//////////////////////////////////////////////////////////////////////////////*/

/*##############################################################################
// Includes.
//############################################################################*/

#include <Libstatmech.h>

/*##############################################################################
// Libstatmech__Fermion__new().
//############################################################################*/

Libstatmech__Fermion
*Libstatmech__Fermion__new(
  const char *s_name,
  double d_rest_mass,
  int i_multiplicity,
  double d_charge
)
{

  Libstatmech__Fermion *self;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if( d_rest_mass < 0. )
    LIBSTATMECH__ERROR( "Rest mass must be >= 0." );

  if( i_multiplicity <= 0 )
    LIBSTATMECH__ERROR( "Multiplicity must be > 0" );

  /*============================================================================
  // Allocate memory for fermion.
  //==========================================================================*/

  self =
    (Libstatmech__Fermion *)
    malloc( sizeof( Libstatmech__Fermion ) );

  /*============================================================================
  // Assign properties.
  //==========================================================================*/

  self->sxName = xmlCharStrdup( s_name );
  self->dRestMass = d_rest_mass;
  self->iMultiplicity = i_multiplicity;
  self->dCharge = d_charge;

  /*============================================================================
  // Initialize standard default quantities.
  //==========================================================================*/

  self->pWorkHash = xmlHashCreate( 0 );

  if( !Libstatmech__value_is_zero( self->dRestMass ) )
  {
    if(
      !Libstatmech__Fermion__updateQuantity(
         self,
         S_NUMBER_DENSITY,
         NULL,
         (Libstatmech__Fermion__Integrand)
            Libstatmech__Fermion__defaultNumberDensityIntegrand
      )
    )
      LIBSTATMECH__ERROR( "Couldn't allocate number density" );
  }
  else
  {
    if(
      !Libstatmech__Fermion__updateQuantity(
         self,
         S_NUMBER_DENSITY,
         (Libstatmech__Fermion__Function)
            Libstatmech__Fermion__defaultNumberDensityFunction,
         NULL
      )
    )
      LIBSTATMECH__ERROR( "Couldn't allocate number density" );
  }

  if(
    !Libstatmech__Fermion__updateQuantity(
       self,
       S_PRESSURE,
       NULL,
       (Libstatmech__Fermion__Integrand)
          Libstatmech__Fermion__defaultPressureIntegrand
    )
  )
    LIBSTATMECH__ERROR( "Couldn't allocate pressure" );

  if(
    !Libstatmech__Fermion__updateQuantity(
       self,
       S_ENERGY_DENSITY,
       NULL,
       (Libstatmech__Fermion__Integrand)
          Libstatmech__Fermion__defaultEnergyDensityIntegrand
    )
  )
    LIBSTATMECH__ERROR( "Couldn't allocate energy density" );

  if(
    !Libstatmech__Fermion__updateQuantity(
       self,
       S_INTERNAL_ENERGY_DENSITY,
       NULL,
       (Libstatmech__Fermion__Integrand)
          Libstatmech__Fermion__defaultInternalEnergyDensityIntegrand
    )
  )
    LIBSTATMECH__ERROR( "Couldn't allocate internal energy density" );

  if(
    !Libstatmech__Fermion__updateQuantity(
       self,
       S_ENTROPY_DENSITY,
       NULL,
       (Libstatmech__Fermion__Integrand)
          Libstatmech__Fermion__defaultEntropyDensityIntegrand
    )
  )
    LIBSTATMECH__ERROR( "Couldn't allocate entropy density" );

  /*============================================================================
  // Done.
  //==========================================================================*/

  return self;

}

/*##############################################################################
// Libstatmech__Fermion__getName().
//############################################################################*/

const char *
Libstatmech__Fermion__getName(
  const Libstatmech__Fermion *self
)
{

  if( !self )
    LIBSTATMECH__ERROR( "Invalid input fermion" );

  return (char *) self->sxName;

}

/*##############################################################################
// Libstatmech__Fermion__getRestMass().
//############################################################################*/

double 
Libstatmech__Fermion__getRestMass(
  const Libstatmech__Fermion *self
)
{

  if( !self )
    LIBSTATMECH__ERROR( "Invalid input fermion" );

  return self->dRestMass;

}

/*##############################################################################
// Libstatmech__Fermion__getMultiplicity().
//############################################################################*/

int 
Libstatmech__Fermion__getMultiplicity(
  const Libstatmech__Fermion *self
)
{

  if( !self )
    LIBSTATMECH__ERROR( "Invalid input fermion" );

  return self->iMultiplicity;

}


/*##############################################################################
// Libstatmech__Fermion__getCharge().
//############################################################################*/

double 
Libstatmech__Fermion__getCharge(
  const Libstatmech__Fermion *self
)
{

  if( !self )
    LIBSTATMECH__ERROR( "Invalid input fermion" );

  return self->dCharge;

}

/*##############################################################################
// Libstatmech__Fermion__free().
//############################################################################*/

void
Libstatmech__Fermion__free(
  Libstatmech__Fermion *self
)
{ 
 
  xmlHashFree(
    self->pWorkHash,
    (xmlHashDeallocator) Libstatmech__Work__free
  );

  xmlFree( self->sxName );
  free( self );

}

/*##############################################################################
// Libstatmech__Fermion__computeIntegrandValue().
//############################################################################*/

double
Libstatmech__Fermion__computeIntegrandValue(
  Libstatmech__Fermion *self,
  const char *s_quantity_name,
  double d_x,
  double d_T,
  double d_mu_kT,
  void *p_user_data
)
{

  Libstatmech__Work *p_work;

  p_work =
    (Libstatmech__Work *)
    xmlHashLookup(
      self->pWorkHash,
      (const xmlChar *) s_quantity_name
    );

  if( !p_work )
    LIBSTATMECH__ERROR( "No such integrand" );

  return
    p_work->pfIntegrand(
      self, d_x, d_T, d_mu_kT, p_user_data
    );

}

/*##############################################################################
// Libstatmech__Fermion__computeQuantity().
//############################################################################*/

double
Libstatmech__Fermion__computeQuantity(
  Libstatmech__Fermion *self,
  const char *s_quantity_name,
  double d_T,
  double d_mu_kT,
  void *p_function_data,
  void *p_integrand_data
)
{

  double d_result;
  Libstatmech__Work *p_work;

  if( !self )
    LIBSTATMECH__ERROR( "Invalid input" );

  if( d_T <= 0. )
    LIBSTATMECH__ERROR( "Temperature must be > 0" );

  p_work =
    (Libstatmech__Work *)
    xmlHashLookup( self->pWorkHash, (const xmlChar *) s_quantity_name );

  if( !p_work ) LIBSTATMECH__ERROR( "No such quantity" );

  p_work->dAlpha = d_mu_kT;
  p_work->dT = d_T;
  p_work->pFermion = self;
  p_work->pBoson = NULL;

  p_work->pFunctionData = p_function_data;
  p_work->pIntegrandData = p_integrand_data;

  d_result =
    Libstatmech__Fermion__integrator( p_work );

  return d_result;

}

/*##############################################################################
// Libstatmech__Fermion__integrator().
//############################################################################*/

double
Libstatmech__Fermion__integrator(
  Libstatmech__Work *p_work
)
{

  double d_result = 0.;

  /*============================================================================
  // First compute value from user-supplied function.
  //==========================================================================*/

  if( p_work->pfFunction )
    d_result =
      p_work->pfFunction(
        p_work->pFermion,
        p_work->dT,
        p_work->dAlpha,
        p_work->pFunctionData
      );

  if( !p_work->pfIntegrand ) return d_result;

  /*============================================================================
  // Now do integral.
  //==========================================================================*/

  if( p_work->dIntegralLowerLimit > p_work->dIntegralUpperLimit )
    LIBSTATMECH__ERROR( "Integral lower limit is larger than upper limit" );

  if(
    p_work->dAlpha > p_work->dIntegralLowerLimit &&
    p_work->dAlpha < p_work->dIntegralUpperLimit
  )
  {
    d_result +=
      Libstatmech__integrate_2(
        p_work,
        p_work->dIntegralLowerLimit,
        p_work->dAlpha
      );
    d_result +=
      Libstatmech__integrate_1(
        p_work,
        p_work->dAlpha,
        p_work->dIntegralUpperLimit
      );
  }
  else
    d_result +=
      Libstatmech__integrate_1(
        p_work,
        p_work->dIntegralLowerLimit,
        p_work->dIntegralUpperLimit
      );

  return d_result;
 
}

/*##############################################################################
// Libstatmech__integrate_1().
//############################################################################*/

double
Libstatmech__integrate_1(
  Libstatmech__Work *p_work,
  double d_start,
  double d_end
)
{

  size_t i_limit = I_LIMIT;
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc( I_WORKSPACE );
       
  double d_result, d_error;

  gsl_function F;
  F.function = &Libstatmech__integrate_helper;
  F.params = p_work;

  /*============================================================================
  // Perform integral.
  //==========================================================================*/

  if( gsl_isinf( p_work->dIntegralUpperLimit ) == 1 )
  {
    gsl_integration_qagiu(
      &F,
      d_start,
      p_work->dEpsAbsolute,
      p_work->dEpsRelative,
      i_limit,
      w,
      &d_result,
      &d_error
    ); 
  }
  else
  {
    gsl_integration_qags(
      &F,
      d_start,
      d_end,
      p_work->dEpsAbsolute,
      p_work->dEpsRelative,
      i_limit,
      w,
      &d_result,
      &d_error
    ); 
  }
     
  gsl_integration_workspace_free( w );

  return d_result; 

}

/*##############################################################################
// Libstatmech__integrate_2().
//############################################################################*/

double
Libstatmech__integrate_2(
  Libstatmech__Work *p_work,
  double d_start,
  double d_end
)
{

  size_t i_limit = I_LIMIT;
  gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc( I_WORKSPACE );
       
  double d_result, d_error;

  gsl_function F;
  F.function = &Libstatmech__integrate_helper;
  F.params = p_work;

  /*============================================================================
  // Perform integral.
  //==========================================================================*/

  gsl_integration_qags(
    &F,
    d_start,
    d_end,
    p_work->dEpsAbsolute,
    p_work->dEpsRelative,
    i_limit,
    w,
    &d_result,
    &d_error
  ); 
     
  gsl_integration_workspace_free( w );

  return d_result; 

}

/*##############################################################################
// Libstatmech__integrate_helper().
//############################################################################*/

double
Libstatmech__integrate_helper(
  double d_x,
  void *p_params
)
{

  Libstatmech__Work *p_work = ( Libstatmech__Work * ) p_params;

  if( p_work->pFermion ) {

    return
      p_work->pfIntegrand(
        p_work->pFermion,
        d_x,
        p_work->dT,
        p_work->dAlpha,
        p_work->pIntegrandData
      );

  }

  if( p_work->pBoson ) {

    return
      p_work->pfIntegrand(
        p_work->pBoson,
        d_x,
        p_work->dT,
        p_work->dAlpha,
        p_work->pIntegrandData
      );

  }

  LIBSTATMECH__ERROR( "Should not reach this point" );

}

/*##############################################################################
// Libstatmech__Fermion__defaultPressureIntegrand().
//############################################################################*/

double
Libstatmech__Fermion__defaultPressureIntegrand(
  Libstatmech__Fermion *self,
  double d_x,
  double d_T,
  double d_alpha,
  void *p_user_data
){

  double d_f, d_mc2, d_gamma, d_part1, d_part2;

  if( p_user_data ) LIBSTATMECH__ERROR( "No user data in default integrand" );

  d_mc2 =
    self->dRestMass * GSL_CONST_CGSM_ELECTRON_VOLT * GSL_CONST_NUM_MEGA;

  d_gamma =
    d_mc2 / ( GSL_CONST_CGSM_BOLTZMANN * d_T );

  /*----------------------------------------------------------------------------
  // Compute part 1.
  //--------------------------------------------------------------------------*/
 
  if( d_alpha - d_x <= 0 )
    d_part1 =
      gsl_log1p( 
        Libstatmech__exp( d_alpha - d_x )
      );
  
  else 
    d_part1 =
      d_alpha - d_x +
      gsl_log1p(
        Libstatmech__exp( d_x - d_alpha )
      );

  /*----------------------------------------------------------------------------
  // Compute part 2.
  //--------------------------------------------------------------------------*/

  if( d_x + 2. * d_gamma + d_alpha <= 0. )
    d_part2 =
      -d_x - 2. * d_gamma - d_alpha +
      gsl_log1p(
        Libstatmech__exp( d_x + 2. * d_gamma + d_alpha )
      );

  else
    d_part2 =
      gsl_log1p(
        Libstatmech__exp( -d_x - 2. * d_gamma - d_alpha )
      );

  d_f =
    sqrt ( d_x * d_x + 2 * d_x * d_gamma ) *
    ( d_x + d_gamma ) *
    ( d_part1 + d_part2 ); 

  d_f *=
    gsl_pow_4( GSL_CONST_CGSM_BOLTZMANN * d_T  ) *
    Libstatmech__Fermion__getMultiplicity( self ) /
    (
      2. *
      gsl_pow_2( M_PI ) *
      gsl_pow_3(
        GSL_CONST_CGSM_PLANCKS_CONSTANT_HBAR *
        GSL_CONST_CGSM_SPEED_OF_LIGHT
      )
    );

  return d_f;

}

/*##############################################################################
// Libstatmech__Fermion__defaultNumberDensityIntegrand().
//############################################################################*/

double
Libstatmech__Fermion__defaultNumberDensityIntegrand(
  Libstatmech__Fermion *self,
  double d_x,
  double d_T,
  double d_alpha,
  void *p_user_data
){

  double d_f, d_gamma, d_mc2, d_part1, d_part2;

  if( p_user_data ) LIBSTATMECH__ERROR( "No user data in default function" );

  d_mc2 =
    self->dRestMass * GSL_CONST_CGSM_ELECTRON_VOLT * GSL_CONST_NUM_MEGA;

  d_gamma =
    d_mc2 / ( GSL_CONST_CGSM_BOLTZMANN * d_T );

  /*----------------------------------------------------------------------------
  // Compute part 1.
  //--------------------------------------------------------------------------*/

  if( d_x - d_alpha >= 0 )
    d_part1 =
      Libstatmech__exp( d_alpha - d_x ) /
      ( Libstatmech__exp( d_alpha - d_x ) + 1. );
  else
    d_part1 =
      1. / ( Libstatmech__exp( d_x - d_alpha ) + 1. );

  /*----------------------------------------------------------------------------
  // Compute part 2.
  //--------------------------------------------------------------------------*/

  if( d_x + 2. * d_gamma + d_alpha >= 0 )
    d_part2 =
      Libstatmech__exp(
        -( d_x + 2. * d_gamma + d_alpha )
      ) /
      (
        Libstatmech__exp(
          -( d_x + 2. * d_gamma + d_alpha )
        ) +
        1.
      );
  else
    d_part2 =
      1. / ( Libstatmech__exp( d_x + 2. * d_gamma + d_alpha ) + 1. );

  /*----------------------------------------------------------------------------
  // Compute integrand.  Check for case that -alpha nearly = gamma.
  //--------------------------------------------------------------------------*/

  if( !gsl_fcmp( d_gamma, -d_alpha, D_GAMMA_ALPHA ) )
    d_f =
      sqrt ( d_x * d_x + 2 * d_x * d_gamma ) *
      ( d_x + d_gamma ) *
      gsl_expm1( 2. * ( d_alpha + d_gamma ) ) * d_part1 * d_part2;
  else
    d_f =
      sqrt ( d_x * d_x + 2 * d_x * d_gamma ) *
      ( d_x + d_gamma ) * ( d_part1 - d_part2 );

  /*----------------------------------------------------------------------------
  // Apply factor in front.
  //--------------------------------------------------------------------------*/

  d_f *=
    gsl_pow_3( GSL_CONST_CGSM_BOLTZMANN * d_T ) *
    Libstatmech__Fermion__getMultiplicity( self ) /
    (
      2. *
      gsl_pow_2( M_PI ) *
      gsl_pow_3(
        GSL_CONST_CGSM_PLANCKS_CONSTANT_HBAR *
        GSL_CONST_CGSM_SPEED_OF_LIGHT
      )
    );

  return d_f;

}

/*##############################################################################
// Libstatmech__Fermion__defaultNumberDensityFunction().
//############################################################################*/

double
Libstatmech__Fermion__defaultNumberDensityFunction(
  Libstatmech__Fermion *self,
  double d_T,
  double d_alpha,
  void *p_user_data
){

  if( p_user_data ) LIBSTATMECH__ERROR( "No user data in default function" );

  return
    (
      gsl_pow_3( M_PI ) * d_alpha / 3.
      +
      gsl_pow_3( d_alpha ) / 3.
    ) *
    gsl_pow_3( GSL_CONST_CGSM_BOLTZMANN * d_T ) *
    Libstatmech__Fermion__getMultiplicity( self ) /
    (
      2. *
      gsl_pow_2( M_PI ) *
      gsl_pow_3(
        GSL_CONST_CGSM_PLANCKS_CONSTANT_HBAR *
        GSL_CONST_CGSM_SPEED_OF_LIGHT
      )
    );

}

/*##############################################################################
// Libstatmech__Fermion__defaultEnergyDensityIntegrand().
//############################################################################*/

double
Libstatmech__Fermion__defaultEnergyDensityIntegrand(
  Libstatmech__Fermion *self,
  double d_x,
  double d_T,
  double d_alpha,
  void *p_user_data
){

  double d_gamma, d_mc2, d_f;

  d_mc2 =
    self->dRestMass * GSL_CONST_CGSM_ELECTRON_VOLT * GSL_CONST_NUM_MEGA;

  d_gamma =
    d_mc2 / ( GSL_CONST_CGSM_BOLTZMANN * d_T );

  d_f =
    GSL_CONST_CGSM_BOLTZMANN * d_T *
    ( d_x + d_gamma ) *
    Libstatmech__Fermion__defaultNumberDensityIntegrand(
      self,
      d_x,
      d_T,
      d_alpha,
      p_user_data
    );

  return d_f;

}

/*##############################################################################
// Libstatmech__Fermion__defaultInternalEnergyDensityIntegrand().
//############################################################################*/

double
Libstatmech__Fermion__defaultInternalEnergyDensityIntegrand(
  Libstatmech__Fermion *self,
  double d_x,
  double d_T,
  double d_alpha,
  void *p_user_data
){

  if( p_user_data ) LIBSTATMECH__ERROR( "No user data in default function" );

  return
    GSL_CONST_CGSM_BOLTZMANN * 
    d_T *
    d_x *
    Libstatmech__Fermion__defaultNumberDensityIntegrand(
      self,
      d_x,
      d_T,
      d_alpha,
      p_user_data
    );

}

/*##############################################################################
// Libstatmech__Fermion__defaultEntropyDensityIntegrand().
//############################################################################*/

double
Libstatmech__Fermion__defaultEntropyDensityIntegrand(
  Libstatmech__Fermion *self,
  double d_x,
  double d_T,
  double d_alpha,
  void *p_user_data
){

  double d_f, d_gamma, d_mc2, d_part1, d_part2, d_part3, d_part4;

  if( p_user_data ) LIBSTATMECH__ERROR( "No user data in default function" );

  d_mc2 =
    self->dRestMass * GSL_CONST_CGSM_ELECTRON_VOLT * GSL_CONST_NUM_MEGA;

  d_gamma =
    d_mc2 / ( GSL_CONST_CGSM_BOLTZMANN * d_T );

  /*----------------------------------------------------------------------------
  // Compute part 1.
  //--------------------------------------------------------------------------*/

  if ( d_x - d_alpha >= 0. )
    d_part1 =
      ( d_x - d_alpha ) * Libstatmech__exp( d_alpha - d_x ) /
      ( Libstatmech__exp( d_alpha - d_x ) + 1. );
  else
    d_part1 =
      ( d_x - d_alpha ) /
      ( Libstatmech__exp( d_x - d_alpha ) + 1. );

  /*----------------------------------------------------------------------------
  // Compute part 2.
  //--------------------------------------------------------------------------*/

  if ( d_x - d_alpha >= 0. )
    d_part2 =
      gsl_log1p(
        Libstatmech__exp( d_alpha - d_x )
      );
  else
    d_part2 =
      d_alpha - d_x +
      gsl_log1p(
        Libstatmech__exp( d_x - d_alpha )
      );
 
  /*----------------------------------------------------------------------------
  // Compute part 3.
  //--------------------------------------------------------------------------*/

  if ( d_x + 2. * d_gamma + d_alpha >= 0 )
    d_part3 =
      gsl_log1p(
        Libstatmech__exp( -( d_x + 2. * d_gamma + d_alpha ) )   
      ); 
  else
    d_part3 =
      -( d_x + 2. * d_gamma + d_alpha ) +
      gsl_log1p(
        Libstatmech__exp( d_x + 2. * d_gamma + d_alpha )
      );

  /*----------------------------------------------------------------------------
  // Compute part 4.
  //--------------------------------------------------------------------------*/
  
  if ( d_x + 2. * d_gamma + d_alpha >= 0 )
    d_part4 =
      ( d_x + 2. * d_gamma + d_alpha ) *
      Libstatmech__exp( - d_x - 2. * d_gamma - d_alpha ) /
      (
        1. +
        Libstatmech__exp( -d_x - 2. * d_gamma - d_alpha )
      );
 
  else
    d_part4 =
      ( d_x + 2. * d_gamma + d_alpha ) /
      (
        1. +
        Libstatmech__exp( d_x + 2. * d_gamma + d_alpha )
      );

  d_f =
    sqrt ( d_x * d_x + 2 * d_x * d_gamma ) *
    ( d_x + d_gamma ) *
    ( d_part1 + d_part2 + d_part3 + d_part4 ); 

  d_f *=
    gsl_pow_3( GSL_CONST_CGSM_BOLTZMANN * d_T ) *
    Libstatmech__Fermion__getMultiplicity( self ) *
    GSL_CONST_CGSM_BOLTZMANN  /
    (
      2. *
      gsl_pow_2( M_PI ) *
      gsl_pow_3(
        GSL_CONST_CGSM_PLANCKS_CONSTANT_HBAR *
        GSL_CONST_CGSM_SPEED_OF_LIGHT
      )
    );

  return d_f;

}

/*##############################################################################
// Libstatmech__Fermion__computeChemicalPotential().
//############################################################################*/

double
Libstatmech__Fermion__computeChemicalPotential(
  Libstatmech__Fermion *self,
  double d_T,
  double d_number_density,
  void *p_function_data,
  void *p_integrand_data
)
{

  typedef struct {
    Libstatmech__Fermion *pFermion;
    double dT;
    double dNumberDensity;
    void *pFunctionData;
    void *pIntegrandData;
  } number_density_root_data;

  int i_status;
  int i_iter = 0, i_max_iter = I_ITER_MAX;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double d_mukT, d_lo = -1., d_hi = 1.;
  gsl_function F;
  number_density_root_data *p_number_density_root_data;

  p_number_density_root_data =
    ( number_density_root_data * ) malloc( sizeof( number_density_root_data ) );

  p_number_density_root_data->pFermion = self;
  p_number_density_root_data->dT = d_T;
  p_number_density_root_data->dNumberDensity = d_number_density;
  p_number_density_root_data->pFunctionData = p_function_data;
  p_number_density_root_data->pIntegrandData = p_integrand_data;

  F.function = &Libstatmech__Fermion__numberDensityRootFinder;
  F.params = p_number_density_root_data;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc( T );

  if(
      !Libstatmech__bracket_root_of_function(
         F, &d_lo, &d_hi
      )
  )
     LIBSTATMECH__ERROR( "Couldn't bracket root" );

  gsl_root_fsolver_set( s, &F, d_lo, d_hi );

  do
    {
       i_iter++;
       i_status = gsl_root_fsolver_iterate( s );
       d_mukT = gsl_root_fsolver_root( s );
       d_lo = gsl_root_fsolver_x_lower( s );
       d_hi = gsl_root_fsolver_x_upper( s );
       i_status =
         gsl_root_test_interval( d_lo, d_hi, D_EPS_ROOT, D_EPS_ROOT );
    }
  while( i_status == GSL_CONTINUE && i_iter < i_max_iter );

  if( i_status == GSL_CONTINUE || i_iter > i_max_iter )
    LIBSTATMECH__ERROR( "Couldn't find number density root" );

  gsl_root_fsolver_free( s );
  free( p_number_density_root_data );

  return d_mukT;

}

/*##############################################################################
// Libstatmech__Fermion__numberDensityRootFinder().
//############################################################################*/

double
Libstatmech__Fermion__numberDensityRootFinder(
  double d_x, void *p_params
)
{

  typedef struct {
    Libstatmech__Fermion *pFermion;
    double dT;
    double dNumberDensity;
    void *pFunctionData;
    void *pIntegrandData;
  } data;

  data *p_data;

  p_data = ( data * ) p_params; 

  return
    Libstatmech__Fermion__computeQuantity(
      p_data->pFermion,
      S_NUMBER_DENSITY,
      p_data->dT,
      d_x,
      p_data->pFunctionData,
      p_data->pIntegrandData
    ) -
    p_data->dNumberDensity;

}

/*##############################################################################
// Libstatmech__bracket_of_root_function()
//############################################################################*/

int
Libstatmech__bracket_root_of_function(
  gsl_function F, double *p_x1, double *p_x2
)
{

  double d_f1, d_f2, d_factor = 1.6, d_xT, d_fT;
  int i_iter = 0, i_max_iter = 1000;

  d_f1 = F.function( *p_x1, F.params );
  d_f2 = F.function( *p_x2, F.params );

  while( i_iter < i_max_iter ) {
    i_iter++;
    if( d_f1 * d_f2 < 0 ) return 1;
    if( fabs( d_f1 ) < fabs( d_f2 ) ) {
      d_xT = *p_x1;
      d_fT = d_f1;
      *p_x1 = *p_x1 + d_factor * ( *p_x1 - *p_x2 );
      *p_x2 = d_xT;
      d_f1 = F.function( *p_x1, F.params );
      d_f2 = d_fT;
    } else {
      d_xT = *p_x2;
      d_fT = d_f2;
      *p_x2 = *p_x2 + d_factor * ( *p_x2 - *p_x1 );
      *p_x1 = d_xT;
      d_f2 = F.function( *p_x2, F.params );
      d_f1 = d_fT;
    }
  }

  return 0;

}

/*##############################################################################
// Libstatmech__Fermion__computeTemperatureDerivative().
//############################################################################*/

double
Libstatmech__Fermion__computeTemperatureDerivative(
  Libstatmech__Fermion *self,
  const char *s_function,
  double d_T,
  double d_number_density,
  void *p_function_data,
  void *p_integrand_data
)
{

  gsl_function F;
  double d_result, d_abserr;

  typedef struct temperature_derivative_data {
    Libstatmech__Fermion *pFermion;
    const char *sFunction;
    double dNumberDensity;
    void *pFunctionData;
    void *pIntegrandData;
  } temperature_derivative_data; 

  temperature_derivative_data *p_temperature_derivative_data;

  p_temperature_derivative_data =
    ( temperature_derivative_data * )
    malloc( sizeof ( temperature_derivative_data ) );

  p_temperature_derivative_data->pFermion = self; 
  p_temperature_derivative_data->sFunction = s_function; 
  p_temperature_derivative_data->dNumberDensity = d_number_density; 
  p_temperature_derivative_data->pFunctionData= p_function_data; 
  p_temperature_derivative_data->pIntegrandData = p_integrand_data; 

  F.function = &Libstatmech__Fermion__differentiate_helper;
  F.params = p_temperature_derivative_data; 
  
  gsl_deriv_central( &F, d_T, D_STEP_FRACTION * d_T, &d_result, &d_abserr );

  free( p_temperature_derivative_data );

  return( d_result );
  
}

/*##############################################################################
// Libstatmech__Fermion__differentiate_helper().
//############################################################################*/

double
Libstatmech__Fermion__differentiate_helper(
  double d_x,
  void *p_params
)
{

  typedef struct data {
    Libstatmech__Fermion *pFermion;
    const char *sFunction;
    double dNumberDensity;
    void *pFunctionData;
    void *pIntegrandData;
  } data; 

  data *p_data;
  double d_mukT;

  p_data = ( data * ) p_params; 

  if( !strcmp( p_data->sFunction, S_CHEMICAL_POTENTIAL ) ) {
    
    return
      Libstatmech__Fermion__computeChemicalPotential(
        p_data->pFermion,
        d_x,
        p_data->dNumberDensity, 
        p_data->pFunctionData,
        p_data->pIntegrandData
      );

  }
  
  else {

    d_mukT =
      Libstatmech__Fermion__computeChemicalPotential(
        p_data->pFermion,
        d_x,
        p_data->dNumberDensity, 
        p_data->pFunctionData,
        p_data->pIntegrandData
      );

    return
      Libstatmech__Fermion__computeQuantity(
        p_data->pFermion,
        p_data->sFunction,
        d_x,
        d_mukT,
        p_data->pFunctionData,
        p_data->pIntegrandData
      );

  }

}

/*##############################################################################
// Libstatmech__Fermion__updateQuantity().
//############################################################################*/

int
Libstatmech__Fermion__updateQuantity(
  Libstatmech__Fermion *self,
  const char *s_quantity_name,
  Libstatmech__Fermion__Function pf_function,
  Libstatmech__Fermion__Integrand pf_integrand
)
{

  Libstatmech__Work *p_work;

  if( !pf_function && !pf_integrand )
    LIBSTATMECH__ERROR( "Quantity must have function or integrand" );

  if(
      !(
         p_work = Libstatmech__Work__new()
      )
  )
    return 0;

  /*============================================================================
  // If setting to default integrand, check that standard quantity.
  //==========================================================================*/

  if(
      ( pf_integrand == DEFAULT_INTEGRAND )
      &&
      !Libstatmech__is_standard_quantity( s_quantity_name )
  )
    return 0;

  /*============================================================================
  // Set function and integrand.
  //==========================================================================*/

  p_work->pfFunction = pf_function;  
  
  if( pf_integrand == DEFAULT_INTEGRAND )
    p_work->pfIntegrand =
      Libstatmech__Fermion__get_default_integrand( s_quantity_name );
  else
    p_work->pfIntegrand = pf_integrand;  

  /*============================================================================
  // Update hash and return.
  //==========================================================================*/

  return
    xmlHashUpdateEntry(
      self->pWorkHash,
      (const xmlChar *) s_quantity_name,
      p_work,
      (xmlHashDeallocator) Libstatmech__Work__free
    ) + 1;

}

/*##############################################################################
// Libstatmech__Work__new().
//############################################################################*/

Libstatmech__Work *
Libstatmech__Work__new( void )
{

  Libstatmech__Work * self =
    ( Libstatmech__Work * ) malloc( sizeof( Libstatmech__Work ) );

  if( !self ) LIBSTATMECH__ERROR( "Couldn't allocate memory" );

  /*============================================================================
  // Initialize parameters.
  //==========================================================================*/

  self->pFermion = NULL;
  self->pBoson = NULL;

  self->pfFunction = NULL;
  self->pFunctionData = NULL;
  self->pfIntegrand = NULL;
  self->pIntegrandData = NULL;

  self->dIntegralLowerLimit = 0.;
  self->dIntegralUpperLimit = GSL_POSINF;

  self->dEpsAbsolute = D_EPS_ABSOLUTE;
  self->dEpsRelative = D_EPS_RELATIVE;

  return self;

}

/*##############################################################################
// Libstatmech__Work__free().
//############################################################################*/

void
Libstatmech__Work__free(
  Libstatmech__Work *self, xmlChar *sx_name
)
{

  if( !sx_name ) LIBSTATMECH__ERROR( "No such quantity" );

  free( self );

}

/*##############################################################################
// Libstatmech__Fermion__updateIntegralLowerLimit().
//############################################################################*/

int
Libstatmech__Fermion__updateIntegralLowerLimit(
  Libstatmech__Fermion *self,
  const char *s_quantity_name,
  double d_lower_limit
)
{

  Libstatmech__Work *p_work;

  if( !self ) LIBSTATMECH__ERROR( "Invalid fermion" );

  if(
    !(
       p_work =
         ( Libstatmech__Work * )
         xmlHashLookup( self->pWorkHash, (const xmlChar *) s_quantity_name )
    )
  )
    LIBSTATMECH__ERROR( "No such quantity" );

  p_work->dIntegralLowerLimit = d_lower_limit;
 
  return 1;

}

/*##############################################################################
// Libstatmech__Fermion__updateIntegralUpperLimit().
//############################################################################*/

int
Libstatmech__Fermion__updateIntegralUpperLimit(
  Libstatmech__Fermion *self,
  const char *s_quantity_name,
  double d_upper_limit
)
{

  Libstatmech__Work *p_work;

  if( !self ) LIBSTATMECH__ERROR( "Invalid fermion" );

  if(
    !(
       p_work =
         ( Libstatmech__Work * )
         xmlHashLookup( self->pWorkHash, (const xmlChar *) s_quantity_name )
    )
  )
    LIBSTATMECH__ERROR( "No such quantity" );

  p_work->dIntegralUpperLimit = d_upper_limit;
 
  return 1;

}

/*##############################################################################
// Libstatmech__Fermion__updateQuantityIntegralAccuracy().
//############################################################################*/

int
Libstatmech__Fermion__updateQuantityIntegralAccuracy(
  Libstatmech__Fermion *self,
  const char *s_quantity_name,
  double d_eps_absolute,
  double d_eps_relative
)
{

  Libstatmech__Work *p_work;

  if(
    !(
       p_work =
         ( Libstatmech__Work * )
         xmlHashLookup( self->pWorkHash, (const xmlChar *) s_quantity_name )
    )
  )
    LIBSTATMECH__ERROR( "No such quantity" );

  p_work->dEpsAbsolute = d_eps_absolute;
  p_work->dEpsRelative = d_eps_relative;
 
  return 1;

}

/*##############################################################################
// Libstatmech__Fermion__get_default_integrand().
//############################################################################*/

Libstatmech__Fermion__Integrand
Libstatmech__Fermion__get_default_integrand( const char *s_quantity_name )
{

  if( !strcmp( s_quantity_name, S_NUMBER_DENSITY ) )
    return
      (Libstatmech__Fermion__Integrand)
      Libstatmech__Fermion__defaultNumberDensityIntegrand;
  else if( !strcmp( s_quantity_name, S_PRESSURE ) )
    return
      (Libstatmech__Fermion__Integrand)
      Libstatmech__Fermion__defaultPressureIntegrand;
  else if( !strcmp( s_quantity_name, S_ENERGY_DENSITY ) )
    return
      (Libstatmech__Fermion__Integrand)
      Libstatmech__Fermion__defaultEnergyDensityIntegrand;
  else if( !strcmp( s_quantity_name, S_INTERNAL_ENERGY_DENSITY ) )
    return
      (Libstatmech__Fermion__Integrand)
      Libstatmech__Fermion__defaultInternalEnergyDensityIntegrand;
  else if( !strcmp( s_quantity_name, S_ENTROPY_DENSITY ) )
    return
      (Libstatmech__Fermion__Integrand)
      Libstatmech__Fermion__defaultEntropyDensityIntegrand;
  else
    return NULL;

}
  
/*##############################################################################
// Routines for Bosons.
//############################################################################*/

/*##############################################################################
// Libstatmech__Boson__new().
//############################################################################*/

Libstatmech__Boson
*Libstatmech__Boson__new(
  const char *s_name,
  double d_rest_mass,
  int i_multiplicity,
  double d_charge
)
{

  Libstatmech__Boson *self;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if( d_rest_mass < 0. )
    LIBSTATMECH__ERROR( "Rest mass must be >= 0." );

  if( i_multiplicity <= 0 )
    LIBSTATMECH__ERROR( "Multiplicity must be > 0" );

  /*============================================================================
  // Allocate memory for boson.
  //==========================================================================*/

  self =
    (Libstatmech__Boson *)
    malloc( sizeof( Libstatmech__Boson ) );

  /*============================================================================
  // Assign properties.
  //==========================================================================*/

  self->sxName = xmlCharStrdup( s_name );
  self->dRestMass = d_rest_mass;
  self->iMultiplicity = i_multiplicity;
  self->dCharge = d_charge;

  /*============================================================================
  // Initialize quantities.
  //==========================================================================*/

  self->pWorkHash = xmlHashCreate( 0 );

  if(
    !Libstatmech__Boson__updateQuantity(
       self,
       S_NUMBER_DENSITY,
       NULL,
       (Libstatmech__Boson__Integrand)
          Libstatmech__Boson__defaultNumberDensityIntegrand
    )
  )
    LIBSTATMECH__ERROR( "Couldn't allocate number density" );

  if(
    !Libstatmech__Boson__updateQuantity(
       self,
       S_PRESSURE,
       NULL,
       (Libstatmech__Boson__Integrand)
          Libstatmech__Boson__defaultPressureIntegrand
    )
  )
    LIBSTATMECH__ERROR( "Couldn't allocate pressure" );

  if(
    !Libstatmech__Boson__updateQuantity(
       self,
       S_ENERGY_DENSITY,
       NULL,
       (Libstatmech__Boson__Integrand)
          Libstatmech__Boson__defaultEnergyDensityIntegrand
    )
  )
    LIBSTATMECH__ERROR( "Couldn't allocate energy density" );

  if(
    !Libstatmech__Boson__updateQuantity(
       self,
       S_INTERNAL_ENERGY_DENSITY,
       NULL,
       (Libstatmech__Boson__Integrand)
          Libstatmech__Boson__defaultInternalEnergyDensityIntegrand
    )
  )
    LIBSTATMECH__ERROR( "Couldn't allocate internal energy density" );

  if(
    !Libstatmech__Boson__updateQuantity(
       self,
       S_ENTROPY_DENSITY,
       NULL,
       (Libstatmech__Boson__Integrand)
          Libstatmech__Boson__defaultEntropyDensityIntegrand
    )
  )
    LIBSTATMECH__ERROR( "Couldn't allocate entropy density" );

  /*============================================================================
  // Done.
  //==========================================================================*/

  return self;

}

/*##############################################################################
// Libstatmech__Boson__updateQuantity().
//############################################################################*/

int
Libstatmech__Boson__updateQuantity(
  Libstatmech__Boson *self,
  const char *s_quantity_name,
  Libstatmech__Boson__Function pf_function,
  Libstatmech__Boson__Integrand pf_integrand
)
{

  Libstatmech__Work *p_work;

  if( !pf_function && !pf_integrand )
    LIBSTATMECH__ERROR( "Quantity must have function or integrand" );

  if(
      !(
         p_work = Libstatmech__Work__new()
      )
  )
    return 0;

  /*============================================================================
  // If setting to default integrand, check that standard quantity.
  //==========================================================================*/

  if(
      ( pf_integrand == DEFAULT_INTEGRAND )
      &&
      !Libstatmech__is_standard_quantity( s_quantity_name )
  )
    return 0;

  /*============================================================================
  // Set function and integrand.
  //==========================================================================*/

  p_work->pfFunction = pf_function;  

  if( pf_integrand == DEFAULT_INTEGRAND )
    p_work->pfIntegrand =
      Libstatmech__Boson__get_default_integrand( s_quantity_name );
  else
    p_work->pfIntegrand = pf_integrand;  

  /*============================================================================
  // Update hash and return.
  //==========================================================================*/

  return
    xmlHashUpdateEntry(
      self->pWorkHash,
      (const xmlChar *) s_quantity_name,
      p_work,
      (xmlHashDeallocator) Libstatmech__Work__free
    ) + 1;

}

/*##############################################################################
// Libstatmech__Boson__getName().
//############################################################################*/

const char *
Libstatmech__Boson__getName(
  const Libstatmech__Boson *self
)
{

  if( !self )
    LIBSTATMECH__ERROR( "Invalid input boson" );

  return (char *) self->sxName;

}

/*##############################################################################
// Libstatmech__Boson__getRestMass().
//############################################################################*/

double 
Libstatmech__Boson__getRestMass(
  const Libstatmech__Boson *self
)
{

  if( !self )
    LIBSTATMECH__ERROR( "Invalid input boson" );

  return self->dRestMass;

}

/*##############################################################################
// Libstatmech__Boson__getMultiplicity().
//############################################################################*/

int 
Libstatmech__Boson__getMultiplicity(
  const Libstatmech__Boson *self
)
{

  if( !self )
    LIBSTATMECH__ERROR( "Invalid input boson" );

  return self->iMultiplicity;

}


/*##############################################################################
// Libstatmech__Boson__getCharge().
//############################################################################*/

double 
Libstatmech__Boson__getCharge(
  const Libstatmech__Boson *self
)
{

  if( !self )
    LIBSTATMECH__ERROR( "Invalid input boson" );

  return self->dCharge;

}

/*##############################################################################
// Libstatmech__Boson__free().
//############################################################################*/

void
Libstatmech__Boson__free(
  Libstatmech__Boson *self
)
{ 
 
  xmlHashFree(
    self->pWorkHash,
    (xmlHashDeallocator) Libstatmech__Work__free
  );

  xmlFree( self->sxName );
  free( self );

}

/*##############################################################################
// Libstatmech__Boson__computeIntegrandValue().
//############################################################################*/

double
Libstatmech__Boson__computeIntegrandValue(
  Libstatmech__Boson *self,
  const char *s_quantity_name,
  double d_x,
  double d_T,
  double d_mu_kT,
  void *p_user_data
)
{

  Libstatmech__Work *p_work;

  p_work =
    (Libstatmech__Work *)
    xmlHashLookup(
      self->pWorkHash,
      (const xmlChar *) s_quantity_name
    );

  if( !p_work )
    LIBSTATMECH__ERROR( "No such integrand" );

  return
    p_work->pfIntegrand(
      self, d_x, d_T, d_mu_kT, p_user_data
    );

}

/*##############################################################################
// Libstatmech__Boson__computeQuantity().
//############################################################################*/

double
Libstatmech__Boson__computeQuantity(
  Libstatmech__Boson *self,
  const char *s_quantity_name,
  double d_T,
  double d_mu_kT,
  void *p_function_data,
  void *p_integrand_data
)
{

  double d_result;
  Libstatmech__Work *p_work;

  if( !self )
    LIBSTATMECH__ERROR( "Invalid input" );

  if( d_T <= 0. )
    LIBSTATMECH__ERROR( "Temperature must be > 0" );

  p_work =
    (Libstatmech__Work *)
    xmlHashLookup( self->pWorkHash, (const xmlChar *) s_quantity_name );

  if( !p_work ) LIBSTATMECH__ERROR( "No such quantity" );

  p_work->dAlpha = d_mu_kT;
  p_work->dT = d_T;
  p_work->pBoson = self;
  p_work->pFermion = NULL;

  p_work->pFunctionData = p_function_data;
  p_work->pIntegrandData = p_integrand_data;

  d_result =
    Libstatmech__Boson__integrator( p_work );

  return d_result;

}

/*##############################################################################
// Libstatmech__Boson__integrator().
//############################################################################*/

double
Libstatmech__Boson__integrator(
  Libstatmech__Work *p_work
)
{
  
  double d_result = 0.;

  if(
    p_work->dAlpha >= 0. &&
    !Libstatmech__value_is_zero( p_work->pBoson->dRestMass )
  )
    LIBSTATMECH__ERROR( "alpha for massive boson should be < 0." );

  /*============================================================================
  // First compute value from user-supplied function.
  //==========================================================================*/

  if( p_work->pfFunction )
    d_result =
      p_work->pfFunction(
        p_work->pBoson,
        p_work->dT,
        p_work->dAlpha,
        p_work->pFunctionData
      );

  if( !p_work->pfIntegrand ) return d_result;

  /*============================================================================
  // Now do integral from integral-lower-limit to infinity.
  //==========================================================================*/

  d_result +=
    Libstatmech__integrate_1(
      p_work,
      p_work->dIntegralLowerLimit,
      p_work->dIntegralUpperLimit
    );

  return d_result;
 
}

/*##############################################################################
// Libstatmech__Boson__defaultPressureIntegrand().
//############################################################################*/

double
Libstatmech__Boson__defaultPressureIntegrand(
  Libstatmech__Boson *self,
  double d_x,
  double d_T,
  double d_alpha,
  void *p_user_data
){

  double d_f, d_mc2, d_gamma; 

  if( p_user_data ) LIBSTATMECH__ERROR( "No user data in default integrand" );

  if(
    GSL_SIGN( self->dRestMass ) == GSL_SIGN( -self->dRestMass ) &&
    GSL_SIGN( d_x ) == GSL_SIGN( -d_x )
  )
    return 0;

  d_mc2 =
    self->dRestMass * GSL_CONST_CGSM_ELECTRON_VOLT * GSL_CONST_NUM_MEGA;

  d_gamma =
    d_mc2 / ( GSL_CONST_CGSM_BOLTZMANN * d_T );

  d_f =
    sqrt ( d_x * d_x + 2 * d_x * d_gamma ) *
    ( d_x + d_gamma ) *
    ( -gsl_log1p
      ( 
        -Libstatmech__exp( d_alpha - d_x )
      )
    );

  d_f *=
    gsl_pow_4(
      d_T *
      GSL_CONST_CGSM_BOLTZMANN
    ) *
    Libstatmech__Boson__getMultiplicity( self ) /
    (
      2. *
      gsl_pow_2( M_PI ) *
      gsl_pow_3(
        GSL_CONST_CGSM_PLANCKS_CONSTANT_HBAR *
        GSL_CONST_CGSM_SPEED_OF_LIGHT
      )
    );

  return d_f;

}

/*##############################################################################
// Libstatmech__Boson__defaultNumberDensityIntegrand().
//############################################################################*/

double
Libstatmech__Boson__defaultNumberDensityIntegrand(
  Libstatmech__Boson *self,
  double d_x,
  double d_T,
  double d_alpha,
  void *p_user_data
){

  double d_f, d_gamma, d_mc2;

  if( p_user_data ) LIBSTATMECH__ERROR( "No user data in default function" );

  if(
    GSL_SIGN( self->dRestMass ) == GSL_SIGN( -self->dRestMass ) &&
    GSL_SIGN( d_x ) == GSL_SIGN( -d_x )
  )
    return 0;

  d_mc2 =
    self->dRestMass * GSL_CONST_CGSM_ELECTRON_VOLT * GSL_CONST_NUM_MEGA;

  d_gamma =
    d_mc2 / ( GSL_CONST_CGSM_BOLTZMANN * d_T );

  d_f =
    sqrt ( d_x * d_x + 2 * d_x * d_gamma ) *
    ( d_x + d_gamma ) /
    gsl_expm1( d_x - d_alpha );

  d_f *=
    gsl_pow_3(
      d_T *
      GSL_CONST_CGSM_BOLTZMANN
    ) *
    Libstatmech__Boson__getMultiplicity( self ) /
    (
      2. *
      gsl_pow_2( M_PI ) *
      gsl_pow_3(
        GSL_CONST_CGSM_PLANCKS_CONSTANT_HBAR *
        GSL_CONST_CGSM_SPEED_OF_LIGHT
      )
    );

  return d_f;

}

/*##############################################################################
// Libstatmech__Boson__defaultEnergyDensityIntegrand().
//############################################################################*/

double
Libstatmech__Boson__defaultEnergyDensityIntegrand(
  Libstatmech__Boson *self,
  double d_x,
  double d_T,
  double d_alpha,
  void *p_user_data
){

  double d_gamma, d_mc2;

  if( p_user_data ) LIBSTATMECH__ERROR( "No user data in default function" );

  if(
    GSL_SIGN( self->dRestMass ) == GSL_SIGN( -self->dRestMass ) &&
    GSL_SIGN( d_x ) == GSL_SIGN( -d_x )
  )
    return 0;

  d_mc2 =
    self->dRestMass * GSL_CONST_CGSM_ELECTRON_VOLT * GSL_CONST_NUM_MEGA;

  d_gamma =
    d_mc2 / ( GSL_CONST_CGSM_BOLTZMANN * d_T );

  return
    GSL_CONST_CGSM_BOLTZMANN *
    ( d_x + d_gamma ) *
    d_T *
    Libstatmech__Boson__defaultNumberDensityIntegrand(
      self,
      d_x,
      d_T,
      d_alpha,
      p_user_data
    );

}

/*##############################################################################
// Libstatmech__Boson__defaultInternalEnergyDensityIntegrand().
//############################################################################*/

double
Libstatmech__Boson__defaultInternalEnergyDensityIntegrand(
  Libstatmech__Boson *self,
  double d_x,
  double d_T,
  double d_alpha,
  void *p_user_data
){

  if( p_user_data ) LIBSTATMECH__ERROR( "No user data in default function" );

  if(
    GSL_SIGN( self->dRestMass ) == GSL_SIGN( -self->dRestMass ) &&
    GSL_SIGN( d_x ) == GSL_SIGN( -d_x )
  )
    return 0;

  return
    GSL_CONST_CGSM_BOLTZMANN *
    d_x *
    d_T *
    Libstatmech__Boson__defaultNumberDensityIntegrand(
      self,
      d_x,
      d_T,
      d_alpha,
      p_user_data
    );

}

/*##############################################################################
// Libstatmech__Boson__defaultEnergyDensitySquareIntegrand().
//############################################################################*/

double
Libstatmech__Boson__defaultEnergyDensitySquareIntegrand(
  Libstatmech__Boson *self,
  double d_x,
  double d_T,
  double d_alpha,
  void *p_user_data
){

  double d_f, d_gamma, d_mc2;

  if( p_user_data ) LIBSTATMECH__ERROR( "No user data in default function" );

  if(
    GSL_SIGN( self->dRestMass ) == GSL_SIGN( -self->dRestMass ) &&
    GSL_SIGN( d_x ) == GSL_SIGN( -d_x )
  )
    return 0;

  d_mc2 =
    self->dRestMass * GSL_CONST_CGSM_ELECTRON_VOLT * GSL_CONST_NUM_MEGA;

  d_gamma =
    d_mc2 / ( GSL_CONST_CGSM_BOLTZMANN * d_T );

  d_f =
    ( d_x * d_x + 2 * d_x * d_gamma ) *
    gsl_pow_4( d_x + d_gamma ) / 
    gsl_expm1( d_x - d_alpha );

  d_f *=
    gsl_pow_8( 
      d_T *
      GSL_CONST_CGSM_BOLTZMANN
    ) *
    Libstatmech__Boson__getMultiplicity( self ) /
    (
      2. *
      gsl_pow_4( M_PI ) *
      gsl_pow_6(
        GSL_CONST_CGSM_PLANCKS_CONSTANT_HBAR *
        GSL_CONST_CGSM_SPEED_OF_LIGHT
      ) 
    );

  return d_f;

}

/*##############################################################################
// Libstatmech__Boson__defaultEntropyDensityIntegrand().
//############################################################################*/

double
Libstatmech__Boson__defaultEntropyDensityIntegrand(
  Libstatmech__Boson *self,
  double d_x,
  double d_T,
  double d_alpha,
  void *p_user_data
){

  double d_f, d_gamma, d_mc2, d_part1, d_part2;

  if( p_user_data ) LIBSTATMECH__ERROR( "No user data in default function" );

  d_mc2 =
    self->dRestMass * GSL_CONST_CGSM_ELECTRON_VOLT * GSL_CONST_NUM_MEGA;

  d_gamma =
    d_mc2 / ( GSL_CONST_CGSM_BOLTZMANN * d_T );

  /*----------------------------------------------------------------------------
  // Compute part 1.
  //--------------------------------------------------------------------------*/

  d_part1 =
    gsl_log1p(
      - Libstatmech__exp( d_alpha - d_x )
    );

  /*----------------------------------------------------------------------------
  // Compute part 2.
  //--------------------------------------------------------------------------*/

  d_part2 =
    ( d_alpha - d_x ) *
    Libstatmech__exp( d_alpha - d_x ) /
    ( 
      1. - 
      Libstatmech__exp( d_alpha - d_x )
    );
 
  d_f =
    sqrt ( d_x * d_x + 2 * d_x * d_gamma ) *
    ( d_x + d_gamma ) *
    ( d_part1 + d_part2 );

  d_f *=
    -1 *
    gsl_pow_3( 
      d_T *
      GSL_CONST_CGSM_BOLTZMANN
    ) *
    GSL_CONST_CGSM_BOLTZMANN  *
    Libstatmech__Boson__getMultiplicity( self ) /
    (
      2. *
      gsl_pow_2( M_PI ) *
      gsl_pow_3(
        GSL_CONST_CGSM_PLANCKS_CONSTANT_HBAR *
        GSL_CONST_CGSM_SPEED_OF_LIGHT
      )
    );

  return d_f;

}

/*##############################################################################
// Libstatmech__Boson__computeChemicalPotential().
//############################################################################*/

double
Libstatmech__Boson__computeChemicalPotential(
  Libstatmech__Boson *self,
  double d_T,
  double d_number_density,
  void *p_function_data,
  void *p_integrand_data
)
{

  typedef struct {
    Libstatmech__Boson *pBoson;
    double dT;
    double dNumberDensity;
    void *pFunctionData;
    void *pIntegrandData;
  } number_density_root_data;

  int i_status;
  int i_iter = 0, i_max_iter = I_ITER_MAX;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double d_mukT, d_lo, d_hi;
  gsl_function F;
  number_density_root_data *p_number_density_root_data;

  p_number_density_root_data =
    ( number_density_root_data * ) malloc( sizeof( number_density_root_data ) );

  p_number_density_root_data->pBoson = self;
  p_number_density_root_data->dT = d_T;
  p_number_density_root_data->dNumberDensity = d_number_density;
  p_number_density_root_data->pFunctionData = p_function_data;
  p_number_density_root_data->pIntegrandData = p_integrand_data;

  F.function = &Libstatmech__Boson__numberDensityRootFinder;
  F.params = p_number_density_root_data;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc( T );

  d_hi = D_BOSON_ALPHA_LIMIT;
  d_lo = d_hi - 10.;

  if(
      !Libstatmech__bracket_root_of_function(
         F, &d_lo, &d_hi
      )
  )
     LIBSTATMECH__ERROR( "Couldn't bracket root" );

  gsl_root_fsolver_set( s, &F, d_lo, d_hi );

  do
    {
       i_iter++;
       i_status = gsl_root_fsolver_iterate( s );
       d_mukT = gsl_root_fsolver_root( s );
       d_lo = gsl_root_fsolver_x_lower( s );
       d_hi = gsl_root_fsolver_x_upper( s );
       i_status =
         gsl_root_test_interval( d_lo, d_hi, 0., D_EPS_ROOT );
    }
  while( i_status == GSL_CONTINUE && i_iter < i_max_iter );

  if( i_status == GSL_CONTINUE || i_iter > i_max_iter )
    LIBSTATMECH__ERROR( "Couldn't find number density root" );

  gsl_root_fsolver_free( s );
  free( p_number_density_root_data );

  return d_mukT;

}

/*##############################################################################
// Libstatmech__Boson__numberDensityRootFinder().
//############################################################################*/

double
Libstatmech__Boson__numberDensityRootFinder(
  double d_x, void *p_params
)
{

  typedef struct {
    Libstatmech__Boson *pBoson;
    double dT;
    double dNumberDensity;
    void *pFunctionData;
    void *pIntegrandData;
  } data;

  data *p_data;

  p_data = ( data * ) p_params; 

  return
    Libstatmech__Boson__computeQuantity(
      p_data->pBoson,
      S_NUMBER_DENSITY,
      p_data->dT,
      d_x,
      p_data->pFunctionData,
      p_data->pIntegrandData
    ) -
    p_data->dNumberDensity;

}

/*##############################################################################
// Libstatmech__Boson__computeTemperatureDerivative().
//############################################################################*/

double
Libstatmech__Boson__computeTemperatureDerivative(
  Libstatmech__Boson *self,
  const char *s_function,
  double d_T,
  double d_number_density,
  void *p_function_data,
  void *p_integrand_data
)
{

  gsl_function F;
  double d_result, d_abserr;

  typedef struct temperature_derivative_data {
    Libstatmech__Boson *pBoson;
    const char *sFunction;
    double dNumberDensity;
    void *pFunctionData;
    void *pIntegrandData;
  } temperature_derivative_data; 

  temperature_derivative_data *p_temperature_derivative_data;

  p_temperature_derivative_data =
    ( temperature_derivative_data * )
    malloc( sizeof ( temperature_derivative_data ) );

  p_temperature_derivative_data->pBoson = self; 
  p_temperature_derivative_data->sFunction = s_function; 
  p_temperature_derivative_data->dNumberDensity = d_number_density; 
  p_temperature_derivative_data->pFunctionData= p_function_data; 
  p_temperature_derivative_data->pIntegrandData = p_integrand_data; 

  F.function = &Libstatmech__Boson__differentiate_helper;
  F.params = p_temperature_derivative_data; 
  
  gsl_deriv_central( &F, d_T, D_STEP_FRACTION * d_T, &d_result, &d_abserr );

  free( p_temperature_derivative_data );

  return( d_result );
  
}

/*##############################################################################
// Libstatmech__Boson__differentiate_helper().
//############################################################################*/

double
Libstatmech__Boson__differentiate_helper(
  double d_x,
  void *p_params
)
{

  typedef struct data {
    Libstatmech__Boson *pBoson;
    const char *sFunction;
    double dNumberDensity;
    void *pFunctionData;
    void *pIntegrandData;
  } data; 

  data *p_data;
  double d_mukT;

  p_data = ( data * ) p_params; 

  if( !strcmp( p_data->sFunction, S_CHEMICAL_POTENTIAL ) ) {
    
    return
      Libstatmech__Boson__computeChemicalPotential(
        p_data->pBoson,
        d_x,
        p_data->dNumberDensity, 
        p_data->pFunctionData,
        p_data->pIntegrandData
      );

  }
  
  else {

    d_mukT =
      Libstatmech__Boson__computeChemicalPotential(
        p_data->pBoson,
        d_x,
        p_data->dNumberDensity, 
        p_data->pFunctionData,
        p_data->pIntegrandData
      );

    return
      Libstatmech__Boson__computeQuantity(
        p_data->pBoson,
        p_data->sFunction,
        d_x,
        d_mukT,
        p_data->pFunctionData,
        p_data->pIntegrandData
      );

  }

}

/*##############################################################################
// Libstatmech__Boson__updateIntegralLowerLimit().
//############################################################################*/

int
Libstatmech__Boson__updateIntegralLowerLimit(
  Libstatmech__Boson *self,
  const char *s_quantity_name,
  double d_lower_limit
)
{

  Libstatmech__Work *p_work;

  if( !self ) LIBSTATMECH__ERROR( "Invalid boson" );

  if(
    !(
       p_work =
         ( Libstatmech__Work * )
         xmlHashLookup( self->pWorkHash, (const xmlChar *) s_quantity_name )
    )
  )
    LIBSTATMECH__ERROR( "No such function" );

  p_work->dIntegralLowerLimit = d_lower_limit;

  return 1;

}

/*##############################################################################
// Libstatmech__Boson__updateIntegralUpperLimit().
//############################################################################*/

int
Libstatmech__Boson__updateIntegralUpperLimit(
  Libstatmech__Boson *self,
  const char *s_quantity_name,
  double d_upper_limit
)
{

  Libstatmech__Work *p_work;

  if( !self ) LIBSTATMECH__ERROR( "Invalid boson" );

  if(
    !(
       p_work =
         ( Libstatmech__Work * )
         xmlHashLookup( self->pWorkHash, (const xmlChar *) s_quantity_name )
    )
  )
    LIBSTATMECH__ERROR( "No such function" );

  p_work->dIntegralUpperLimit = d_upper_limit;

  return 1;

}

/*##############################################################################
// Libstatmech__exp().
//############################################################################*/

double
Libstatmech__exp( double d_x )
{

  if( d_x < -700. )
    return 0.;
  else if( d_x > 700. )
    LIBSTATMECH__ERROR( "Too large an argument" );
  else if( fabs( d_x ) < 1.e-10 )
    return ( 1 + d_x );
  else
    return exp( d_x );

}

/*##############################################################################
// Libstatmech__is_standard_quantity().
//############################################################################*/

int
Libstatmech__is_standard_quantity( const char *s_quantity_name )
{

  if( !strcmp( s_quantity_name, S_NUMBER_DENSITY ) )
    return 1;
  else if( !strcmp( s_quantity_name, S_PRESSURE ) )
    return 1;
  else if( !strcmp( s_quantity_name, S_ENERGY_DENSITY ) )
    return 1;
  else if( !strcmp( s_quantity_name, S_INTERNAL_ENERGY_DENSITY ) )
    return 1;
  else if( !strcmp( s_quantity_name, S_ENTROPY_DENSITY ) )
    return 1;
  else
    return 0;

}
  
/*##############################################################################
// Libstatmech__Boson__updateQuantityIntegralAccuracy().
//############################################################################*/

int
Libstatmech__Boson__updateQuantityIntegralAccuracy(
  Libstatmech__Boson *self,
  const char *s_quantity_name,
  double d_eps_absolute,
  double d_eps_relative
)
{

  Libstatmech__Work *p_work;

  if(
    !(
       p_work =
         ( Libstatmech__Work * )
         xmlHashLookup( self->pWorkHash, (const xmlChar *) s_quantity_name )
    )
  )
    LIBSTATMECH__ERROR( "No such function" );

  p_work->dEpsAbsolute = d_eps_absolute;
  p_work->dEpsRelative = d_eps_relative;
 
  return 1;

}

/*##############################################################################
// Libstatmech__Boson__get_default_integrand().
//############################################################################*/

Libstatmech__Boson__Integrand
Libstatmech__Boson__get_default_integrand( const char *s_quantity_name )
{

  if( !strcmp( s_quantity_name, S_NUMBER_DENSITY ) )
    return
      (Libstatmech__Boson__Integrand)
      Libstatmech__Boson__defaultNumberDensityIntegrand;
  else if( !strcmp( s_quantity_name, S_PRESSURE ) )
    return
      (Libstatmech__Boson__Integrand)
      Libstatmech__Boson__defaultPressureIntegrand;
  else if( !strcmp( s_quantity_name, S_ENERGY_DENSITY ) )
    return
      (Libstatmech__Boson__Integrand)
      Libstatmech__Boson__defaultEnergyDensityIntegrand;
  else if( !strcmp( s_quantity_name, S_INTERNAL_ENERGY_DENSITY ) )
    return
      (Libstatmech__Boson__Integrand)
      Libstatmech__Boson__defaultInternalEnergyDensityIntegrand;
  else if( !strcmp( s_quantity_name, S_ENTROPY_DENSITY ) )
    return
      (Libstatmech__Boson__Integrand)
      Libstatmech__Boson__defaultEntropyDensityIntegrand;
  else
    return NULL;

}

/*##############################################################################
// DEFAULT_INTEGRAND(). 
//############################################################################*/

double
DEFAULT_INTEGRAND( void *p_user, double x, double y, double z, void *p_data )
{

  if(
    !p_user ||
    Libstatmech__value_is_zero( x ) ||
    Libstatmech__value_is_zero( y ) ||
    Libstatmech__value_is_zero( z ) ||
    p_data
  )
    LIBSTATMECH__ERROR( "Invalid input" );

  LIBSTATMECH__ERROR( "Dummy function" );

}

/*##############################################################################
// Libstatmech__value_is_zero().
//############################################################################*/

int
Libstatmech__value_is_zero( double d_x )
{

  if( GSL_SIGN( d_x ) == GSL_SIGN( -d_x ) )
    return 1;
  else
    return 0;

}

/*##############################################################################
// Libstatmech__Fermion__copy().
//############################################################################*/

Libstatmech__Fermion *
Libstatmech__Fermion__copy(
  Libstatmech__Fermion * self
)
{

  Libstatmech__Fermion * p_new_fermion;

  if( !self ) LIBSTATMECH__ERROR( "Invalid fermion" );

  p_new_fermion =
    Libstatmech__Fermion__new(
      Libstatmech__Fermion__getName( self ), 
      Libstatmech__Fermion__getRestMass( self ), 
      Libstatmech__Fermion__getMultiplicity( self ), 
      Libstatmech__Fermion__getCharge( self )
    );

  return p_new_fermion;

}

/*##############################################################################
// Libstatmech__Boson__copy().
//############################################################################*/

Libstatmech__Boson *
Libstatmech__Boson__copy(
  Libstatmech__Boson * self
)
{

  Libstatmech__Boson * p_new_boson;

  if( !self ) LIBSTATMECH__ERROR( "Invalid boson" );

  p_new_boson =
    Libstatmech__Boson__new(
      Libstatmech__Boson__getName( self ), 
      Libstatmech__Boson__getRestMass( self ), 
      Libstatmech__Boson__getMultiplicity( self ), 
      Libstatmech__Boson__getCharge( self )
    );

  return p_new_boson;

}
