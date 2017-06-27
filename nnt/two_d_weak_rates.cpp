////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2012 Clemson University.
//
// This file was originally written by Bradley S. Meyer.
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this software; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
// USA
//
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//!
//! \file
//! \brief A file containing routines to compute weak rates depending
//!        on temperature and mass density of electrons.
//!
////////////////////////////////////////////////////////////////////////////////

#include "nnt/two_d_weak_rates.h"

namespace nnt
{

//##############################################################################
// ffnIV_I()
//##############################################################################

double
ffnIV_I(
  double d_eta_F,
  double d_eta_L,
  double d_zeta_n,
  double d_eta_nu
)
{

  if( d_eta_nu == GSL_NEGINF )
  {
    return
      ffnIV__fermi_dirac( 4, d_eta_F - d_eta_L ) 
      +
      ( 2. * d_zeta_n + 4. * d_eta_L ) *
      ffnIV__fermi_dirac( 3, d_eta_F - d_eta_L ) 
      +
      ( 
	6. * gsl_pow_2( d_eta_L ) + 
	6. * d_eta_L * d_zeta_n + 
	gsl_pow_2( d_zeta_n ) 
      ) *
      ffnIV__fermi_dirac( 2, d_eta_F - d_eta_L ) 
      +
      (
	4. * gsl_pow_3( d_eta_L ) +
	6. * gsl_pow_2( d_eta_L ) * d_zeta_n +
	2. * d_eta_L * gsl_pow_2( d_zeta_n )
      ) *
      ffnIV__fermi_dirac( 1, d_eta_F - d_eta_L ) 
      +
      ( 
	gsl_pow_4( d_eta_L ) +
	2. * d_zeta_n * gsl_pow_3( d_eta_L ) +
	gsl_pow_2( d_zeta_n ) * gsl_pow_2( d_eta_L )
      ) * 
      ffnIV__fermi_dirac( 0, d_eta_F - d_eta_L ); 
  }
  else
  {
    return
      ( 
	ffnIV__fermi_dirac( 4, d_eta_F - d_eta_L ) -
	ffnIV__fermi_dirac( 4, d_eta_nu - d_zeta_n - d_eta_L )
      )
      +
      ( 2. * d_zeta_n + 4. * d_eta_L ) *
      ( 
	ffnIV__fermi_dirac( 3, d_eta_F - d_eta_L ) -
	ffnIV__fermi_dirac( 3, d_eta_nu - d_zeta_n - d_eta_L )
      )
      +
      ( 
	6. * gsl_pow_2( d_eta_L ) + 
	6. * d_eta_L * d_zeta_n + 
	gsl_pow_2( d_zeta_n ) 
      ) *
      ( 
	ffnIV__fermi_dirac( 2, d_eta_F - d_eta_L ) -
	ffnIV__fermi_dirac( 2, d_eta_nu - d_zeta_n - d_eta_L )
      )
      +
      (
	4. * gsl_pow_3( d_eta_L ) +
	6. * gsl_pow_2( d_eta_L ) * d_zeta_n +
	2. * d_eta_L * gsl_pow_2( d_zeta_n )
      ) *
      ( 
	ffnIV__fermi_dirac( 1, d_eta_F - d_eta_L ) -
	ffnIV__fermi_dirac( 1, d_eta_nu - d_zeta_n - d_eta_L )
      )
      +
      ( 
	gsl_pow_4( d_eta_L ) +
	2. * d_zeta_n * gsl_pow_3( d_eta_L ) +
	gsl_pow_2( d_zeta_n ) * gsl_pow_2( d_eta_L )
      ) * 
      ( 
	ffnIV__fermi_dirac( 0, d_eta_F - d_eta_L ) -
	ffnIV__fermi_dirac( 0, d_eta_nu - d_zeta_n - d_eta_L )
      );
  }

}
    
//##############################################################################
// ffnIV_Ie()
//##############################################################################

double
ffnIV_Ie(
  double d_mekT,
  double d_eta_F,
  double d_eta_L,
  double d_zeta_n,
  double d_eta_nu
)
{

  if( d_eta_nu == GSL_NEGINF )
    return
      ffnIV_I( d_eta_F, d_eta_L, d_zeta_n, d_eta_nu ) /
      gsl_pow_5( d_mekT );
  else
    return
      ffnIV_I( d_eta_F, d_eta_L, d_zeta_n, d_eta_nu ) /
      ( 
	gsl_pow_5( d_mekT ) *
	( 1. - exp( d_eta_nu - d_zeta_n - d_eta_F ) )
      );

}

//##############################################################################
// ffnIV__fermi_dirac().
//##############################################################################

double
ffnIV__fermi_dirac( int i_n, double d_eta )
{

  double d_n = (double) i_n;

  switch( i_n ) {

    case 0:
      if( d_eta > 500. )
	return d_eta;
      else
	return
	  log( 1 + exp( d_eta ) );

    case 1:
      if( d_eta <= 0 )
	return
	  exp( d_eta );
      else
	return
	  pow( d_eta, d_n + 1. ) / ( d_n + 1. ) + 2. - exp( -d_eta ); 

    case 2:
      if( d_eta <= 0 )
	return
	  gsl_sf_gamma( d_n + 1. ) * exp( d_eta );
      else
	return
	  pow( d_eta, d_n + 1. ) / ( d_n + 1. ) + 
	  2. * gsl_sf_gamma( d_n + 1. ) * d_eta + 
	  gsl_sf_gamma( d_n + 1. ) * exp( -d_eta );

    case 3:
      if( d_eta <= 0 )
	return
	  gsl_sf_gamma( d_n + 1. ) * exp( d_eta );
      else
	return
	  pow( d_eta, d_n + 1. ) / ( d_n + 1. ) + 
	  gsl_pow_2( M_PI ) * gsl_pow_2( d_eta ) / 2. +
	  2. * gsl_sf_gamma( d_n + 1. ) - 
	  gsl_sf_gamma( d_n + 1. ) * exp( -d_eta ); 

    case 4:
      if( d_eta <= 0 )
	return
	  gsl_sf_gamma( d_n + 1. ) * exp( d_eta );
      else
	return
	  pow( d_eta, d_n + 1. ) / ( d_n + 1. ) + 
	  2. * gsl_pow_2( M_PI ) * gsl_pow_3( d_eta ) / 3. +
	  2. * gsl_sf_gamma( d_n + 1. ) * d_eta + 
	  gsl_sf_gamma( d_n + 1. ) * exp( -d_eta ); 

    case 5:
      if( d_eta <= 0 )
	return
	  gsl_sf_gamma( d_n + 1. ) * exp( d_eta );
      else
	return
	  pow( d_eta, d_n + 1. ) / ( d_n + 1. ) + 
	  5. * gsl_pow_2( M_PI ) * gsl_pow_4( d_eta ) / 6. +
	  7. * gsl_pow_2( M_PI ) * gsl_pow_2( d_eta ) / 6. +
	  2. * gsl_sf_gamma( d_n + 1. ) - 
	  gsl_sf_gamma( d_n + 1. ) * exp( -d_eta ); 

    default:
    {
      fprintf(
	stderr,
	"Here the fermi-dirac integral order should be from 0 to 6!\n"
      );
      exit( EXIT_FAILURE );
    }

  }

}

//##############################################################################
// ffnIV__compute_Ie()
//##############################################################################

double
ffnIV__compute_Ie(
  Libnucnet__Reaction *p_reaction,
  Libnucnet__Net *p_net,
  double d_electron_mass,
  double d_t9,
  double d_eta_F,
  double d_mu_nue_kT
)
{

  double d_kT, d_mekT, d_zeta_n, d_eta_L;

  d_kT = compute_kT_in_MeV( d_t9 );

  d_mekT =
    d_electron_mass / d_kT;

  d_zeta_n =
    (
      compute_reaction_nuclear_Qvalue(
	p_net, 
	p_reaction,
	d_electron_mass
      )
    ) / d_kT;

  if( d_zeta_n + d_mekT > 0 )
    d_eta_L = d_mekT;
  else
    d_eta_L = fabs( d_zeta_n );

  if( is_positron_capture_reaction( p_reaction ) )
    return 
      ffnIV_Ie( d_mekT, -d_eta_F, d_eta_L, d_zeta_n, d_mu_nue_kT );
  else
    return 
      ffnIV_Ie( d_mekT, d_eta_F, d_eta_L, d_zeta_n, d_mu_nue_kT );

}

//##############################################################################
// TwoDWeakQuantity::computeValue().
//##############################################################################

std::pair<double,double>
TwoDWeakQuantity::computeValue(
  double d_t9,
  double d_rhoe
)
{

  if( d_t9 <= 0. || d_rhoe < 0. )
  {
    std::cerr << "Invalid input." << std::endl;
    exit( EXIT_FAILURE );
  }

  return
    two_d_interpolation(
      this->pT9Vector,
      this->pLog10RhoeVector,
      this->pMatrix,
      d_t9,
      log10( d_rhoe )
    );

}

//##############################################################################
// TwoDWeakQuantity().
//##############################################################################

TwoDWeakQuantity::TwoDWeakQuantity(
  Libnucnet__Reaction *p_reaction,
  const char *s_property
)
{

  size_t i_rows = 0, i_columns = 0;

  pReaction = p_reaction;

  sReaction = Libnucnet__Reaction__getString( p_reaction );

  pT9Vector = get_rate_function_property_gsl_vector( p_reaction, s_T9 );

  pLog10RhoeVector =
    get_rate_function_property_gsl_vector( p_reaction, s_LOG10_RHOE );

  i_rows = WnMatrix__get_gsl_vector_size( pT9Vector );

  i_columns = WnMatrix__get_gsl_vector_size( pLog10RhoeVector );

  pMatrix =
    gsl_matrix_alloc(
      i_rows,
      i_columns
    );

  Libnucnet__Reaction__iterateUserRateFunctionProperties(
    p_reaction,
    s_property, 
    NULL,
    NULL,
    (Libnucnet__Reaction__user_rate_property_iterate_function)
       get_property_matrix,
    pMatrix
  );

}

//##############################################################################
// ~TwoDWeakQuantity().
//##############################################################################

TwoDWeakQuantity::~TwoDWeakQuantity()
{

  gsl_vector_free( pT9Vector );
  gsl_vector_free( pLog10RhoeVector );
  gsl_matrix_free( pMatrix );

} 

//##############################################################################
// TwoDWeakQuantity() copy.
//##############################################################################

TwoDWeakQuantity::TwoDWeakQuantity( const TwoDWeakQuantity& other )
{

  pReaction = other.pReaction;

  sReaction = other.sReaction;

  pT9Vector = gsl_vector_alloc( other.pT9Vector->size );
  pLog10RhoeVector = gsl_vector_alloc( other.pLog10RhoeVector->size );
  pMatrix =
    gsl_matrix_alloc(
      other.pMatrix->size1,
      other.pMatrix->size2
    );

  gsl_vector_memcpy( pT9Vector, other.pT9Vector );
  gsl_vector_memcpy( pLog10RhoeVector, other.pLog10RhoeVector );
  gsl_matrix_memcpy( pMatrix, other.pMatrix );

}

} // namespace nnt
