//////////////////////////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//!
//! \file
//! \brief Code for user-defined neutrino rate functions.
//!
////////////////////////////////////////////////////////////////////////////////

#include "neutrino_rate_functions.h"

namespace user
{

//##############################################################################
// Global map to store neutrino arrays.
//##############################################################################

boost::unordered_map<std::string, NeutrinoQuantity> my_xsecs;

//##############################################################################
// NeutrinoQuantity::computeValue().
//##############################################################################

double
NeutrinoQuantity::computeValue(
  double d_x
)
{

  return
    nnt::linear_interpolation(
      this->pTVector,
      this->pLog10XSecVector,
      d_x
    );

}

//##############################################################################
// NeutrinoQuantity().
//##############################################################################

NeutrinoQuantity::NeutrinoQuantity(
  Libnucnet__Reaction *p_reaction
)
{

  pReaction = p_reaction;

  sReaction = Libnucnet__Reaction__getString( p_reaction );

  pTVector = nnt::get_rate_function_property_gsl_vector( p_reaction, S_T_NU );

  pLog10XSecVector =
    nnt::get_rate_function_property_gsl_vector( p_reaction, S_LOG10_XSEC );

  if( pTVector->size != pLog10XSecVector->size )
  {
    std::cerr << "Tnu and rate arrays not equal length." << std::endl;
    exit( EXIT_FAILURE );
  }

}

//##############################################################################
// ~NeutrinoQuantity().
//##############################################################################

NeutrinoQuantity::~NeutrinoQuantity()
{

  gsl_vector_free( pTVector );
  gsl_vector_free( pLog10XSecVector );

} 

//##############################################################################
// NeutrinoQuantity() copy.
//##############################################################################

NeutrinoQuantity::NeutrinoQuantity( const NeutrinoQuantity& other )
{

  pReaction = other.pReaction;

  sReaction = other.sReaction;

  pTVector = gsl_vector_alloc( other.pTVector->size );
  pLog10XSecVector = gsl_vector_alloc( other.pLog10XSecVector->size );

  gsl_vector_memcpy( pTVector, other.pTVector );
  gsl_vector_memcpy( pLog10XSecVector, other.pLog10XSecVector );

}

//##############################################################################
// register_my_neutrino_rate_functions().
//##############################################################################

void
register_my_neutrino_rate_functions( Libnucnet__Reac *p_reac )
{

  //============================================================================
  // Register neutrino capture on neutron.
  //============================================================================

  Libnucnet__Reac__registerUserRateFunction(
    p_reac,
    NU_N_CAPTURE,
    (Libnucnet__Reaction__userRateFunction) nu_n_capture_function
  );

  //============================================================================
  // Register anti-neutrino capture on proton.
  //============================================================================

  Libnucnet__Reac__registerUserRateFunction(
    p_reac,
    NU_P_CAPTURE,
    (Libnucnet__Reaction__userRateFunction) nu_p_capture_function
  );

  //============================================================================
  // Register neutrino-nucleus interaction.
  //============================================================================

  Libnucnet__Reac__registerUserRateFunction(
    p_reac,
    NU_NUCL,
    (Libnucnet__Reaction__userRateFunction) nu_nucl_function
  );

}

//##############################################################################
// update_my_neutrino_rate_functions_data().
//##############################################################################

void
update_my_neutrino_rate_functions_data(
  nnt::Zone &zone
)
{

  double d_tau_lum;

  //============================================================================
  // If neutrino luminosity declines with time, first set initial luminosity
  // and then update.
  //============================================================================

  if( zone.hasProperty( nnt::s_TAU_LUM_NEUTRINO ) )
  {
      
    if( zone.hasProperty( S_L_NU_E ) && !zone.hasProperty( S_L_NU_E_0 ) )
      zone.updateProperty(
        S_L_NU_E_0,
        zone.getProperty( S_L_NU_E )
      );

    if(
      zone.hasProperty( S_L_NU_E_BAR ) &&
      !zone.hasProperty( S_L_NU_E_BAR_0 )
    )
      zone.updateProperty(
        S_L_NU_E_BAR_0,
        zone.getProperty( S_L_NU_E_BAR )
      );

    if( zone.hasProperty( S_L_NU_MU ) && !zone.hasProperty( S_L_NU_MU_0 ) )
      zone.updateProperty(
        S_L_NU_MU_0,
        zone.getProperty( S_L_NU_MU )
      );

    if(
      zone.hasProperty( S_L_NU_MU_BAR ) &&
      !zone.hasProperty( S_L_NU_MU_BAR_0 )
    )
      zone.updateProperty(
        S_L_NU_MU_BAR_0,
        zone.getProperty( S_L_NU_MU_BAR )
      );

    if( zone.hasProperty( S_L_NU_TAU ) && !zone.hasProperty( S_L_NU_TAU_0 ) )
      zone.updateProperty(
        S_L_NU_TAU_0,
        zone.getProperty( S_L_NU_TAU )
      );

    if(
      zone.hasProperty( S_L_NU_TAU_BAR ) &&
      !zone.hasProperty( S_L_NU_TAU_BAR_0 )
    )
      zone.updateProperty(
        S_L_NU_TAU_BAR_0,
        zone.getProperty( S_L_NU_TAU_BAR )
      );

    try
    {
      d_tau_lum =
        boost::lexical_cast<double>(
          zone.getProperty( nnt::s_TAU_LUM_NEUTRINO )
        );
    }
    catch( const boost::bad_lexical_cast& e )
    {
      if( zone.getProperty( nnt::s_TAU_LUM_NEUTRINO ) == "inf" )
      {
        d_tau_lum = GSL_POSINF;
      }
      else
      {
	std::cerr << "Invalid tau for neutrino luminosity." << std::endl;
        exit( EXIT_FAILURE );
      }
    }

    double d_factor =
      exp(
        -boost::lexical_cast<double>( zone.getProperty( nnt::s_TIME ) ) /
        d_tau_lum
      );

    zone.updateProperty(
      S_L_NU_E,
      boost::lexical_cast<std::string>(
        boost::lexical_cast<double>( zone.getProperty( S_L_NU_E_0 ) ) *
        d_factor
      )
    );

    zone.updateProperty(
      S_L_NU_E_BAR,
      boost::lexical_cast<std::string>(
        boost::lexical_cast<double>( zone.getProperty( S_L_NU_E_BAR_0 ) ) *
        d_factor
      )
    );

    zone.updateProperty(
      S_L_NU_MU,
      boost::lexical_cast<std::string>(
        boost::lexical_cast<double>( zone.getProperty( S_L_NU_MU_0 ) ) *
        d_factor
      )
    );

    zone.updateProperty(
      S_L_NU_MU_BAR,
      boost::lexical_cast<std::string>(
        boost::lexical_cast<double>( zone.getProperty( S_L_NU_MU_BAR_0 ) ) *
        d_factor
      )
    );

    zone.updateProperty(
      S_L_NU_TAU,
      boost::lexical_cast<std::string>(
        boost::lexical_cast<double>( zone.getProperty( S_L_NU_TAU_0 ) ) *
        d_factor
      )
    );

    zone.updateProperty(
      S_L_NU_TAU_BAR,
      boost::lexical_cast<std::string>(
        boost::lexical_cast<double>( zone.getProperty( S_L_NU_TAU_BAR_0 ) ) *
        d_factor
      )
    );

  }

  //============================================================================
  // Set data for neutrino-nucleus interaction.
  //============================================================================

  if(
    Libnucnet__Reac__isRegisteredRateFunction(
      Libnucnet__Net__getReac(
        Libnucnet__Zone__getNet( zone.getNucnetZone() )
      ),
      NU_NUCL
    )
  )
  {
    Libnucnet__Zone__updateDataForUserRateFunction(
      zone.getNucnetZone(),
      NU_NUCL,
      &zone
    );
  }

  //============================================================================
  // Set data for anti-neutrino capture on proton.
  //============================================================================

  if(
    Libnucnet__Reac__isRegisteredRateFunction(
      Libnucnet__Net__getReac(
        Libnucnet__Zone__getNet( zone.getNucnetZone() )
      ),
      NU_P_CAPTURE
    )
  )
  {
    Libnucnet__Zone__updateDataForUserRateFunction(
      zone.getNucnetZone(),
      NU_P_CAPTURE,
      &zone
    );
  }

  //============================================================================
  // Set data for neutrino capture on neutron.
  //============================================================================

  if(
    Libnucnet__Reac__isRegisteredRateFunction(
      Libnucnet__Net__getReac(
        Libnucnet__Zone__getNet( zone.getNucnetZone() )
      ),
      NU_N_CAPTURE
    )
  )
  {
    Libnucnet__Zone__updateDataForUserRateFunction(
      zone.getNucnetZone(),
      NU_N_CAPTURE,
      &zone
    );
  }

}

//##############################################################################
// set_nu_nucl_hash().
//##############################################################################

void
set_nu_nucl_hash( Libnucnet__Reac * p_reac )
{

  size_t i_size;

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list( p_reac );

  for(
    nnt::reaction_list_t::iterator it = reaction_list.begin();
    it != reaction_list.end();
    it++
  )
  {

    if(
      strcmp(
        Libnucnet__Reaction__getRateFunctionKey( it->getNucnetReaction() ),
        S_NU_NUCL
      ) == 0
    )
    {

      i_size = 0;

      Libnucnet__Reaction__iterateUserRateFunctionProperties(
        it->getNucnetReaction(),
        S_LOG10_XSEC,
        NULL,
        NULL,
        (Libnucnet__Reaction__user_rate_property_iterate_function)
           nnt::get_property_array_size,
        &i_size
      );

      if( i_size > 0 )
      {
        my_xsecs.insert(
          std::make_pair(
            Libnucnet__Reaction__getString( it->getNucnetReaction() ),
            NeutrinoQuantity(
               it->getNucnetReaction()
            )
          )
        );

      }

    }

  }

}
    
//##############################################################################
// nu_nucl_function()
//##############################################################################

double
nu_nucl_function(
  Libnucnet__Reaction *p_reaction,
  double d_t9,
  nnt::Zone& zone
)
{

  char s_neutrino[nnt::i_BUF_SIZE];

  if( d_t9 <= 0 )
  {
    fprintf( stderr, "Must have t9 > 0 for this function.\n" );
    exit( EXIT_FAILURE );
  }

  //============================================================================
  // Get reactant neutrino.
  //============================================================================

  Libnucnet__Reaction__iterateReactants(
    p_reaction,
    (Libnucnet__Reaction__Element__iterateFunction)
       get_reactant_neutrino,
    s_neutrino
  );

  if( strcmp( s_neutrino, "neutrino_e" ) == 0 )
  {
    return
      compute_neutrino_rate(
        zone,
        p_reaction,
        S_T_NU_E,
        S_L_NU_E
      );
  } else if( strcmp( s_neutrino, "anti-neutrino_e" ) == 0 )
  {
    return
      compute_neutrino_rate(
        zone,
        p_reaction,
        S_T_NU_E_BAR,
        S_L_NU_E_BAR
      );
  } else if( strcmp( s_neutrino, "neutrino_mu" ) == 0 )
  {
    return
      compute_neutrino_rate(
        zone,
        p_reaction,
        S_T_NU_MU,
        S_L_NU_MU
      );
  } else if( strcmp( s_neutrino, "anti-neutrino_mu" ) == 0 )
  {
    return
      compute_neutrino_rate(
        zone,
        p_reaction,
        S_T_NU_MU_BAR,
        S_L_NU_MU_BAR
      );
  } else if( strcmp( s_neutrino, "neutrino_tau" ) == 0 )
  {
    return
      compute_neutrino_rate(
        zone,
        p_reaction,
        S_T_NU_TAU,
        S_L_NU_TAU
      );
  } else if( strcmp( s_neutrino, "anti-neutrino_tau" ) == 0 )
  {
    return
      compute_neutrino_rate(
        zone,
        p_reaction,
        S_T_NU_TAU_BAR,
        S_L_NU_TAU_BAR
      );
  } else
  {
    fprintf( stderr, "No such neutrino.\n" );
    exit( EXIT_FAILURE );
  } 

}
    
//##############################################################################
// nu_p_capture_function()
//##############################################################################

double
nu_p_capture_function(
  Libnucnet__Reaction *p_reaction,
  double d_t9,
  nnt::Zone& zone
)
{

  double d_Lnu, d_Enu_MeV, d_prefactor, d_Delta, d_R6; 
  const char *s_tmp;

  if( d_t9  <= 0 )
  {
    fprintf( stderr, "Must have t9 > 0 for this function.\n" );
    exit( EXIT_FAILURE );
  }

  //============================================================================
  // Neutrino parameters.
  //============================================================================

  d_Lnu = boost::lexical_cast<double>( zone.getProperty( S_L_NU_E_BAR ) );

  d_Enu_MeV =
    3.15 * boost::lexical_cast<double>( zone.getProperty( S_T_NU_E_BAR ) );

  d_R6 =
    boost::lexical_cast<double>( zone.getProperty( nnt::s_RADIUS ) ) / 1.e6;

  //============================================================================
  // Rate parameters.
  //============================================================================

  if(
    ( s_tmp =
        Libnucnet__Reaction__getUserRateFunctionProperty(
          p_reaction,
          S_PREFACTOR,
          NULL,
          NULL
        )
    )
  )
    d_prefactor = atof( s_tmp );
  else
  {
    fprintf( stderr, " prefactor not provided.\n" );
    exit( EXIT_FAILURE );
  }

  if(
    ( s_tmp =
        Libnucnet__Reaction__getUserRateFunctionProperty(
          p_reaction,
          S_DELTA,
          NULL,
          NULL
        )
    )
  )
    d_Delta = atof( s_tmp );
  else
  {
    fprintf( stderr, " Delta not provided.\n" );
    exit( EXIT_FAILURE );
  }

  //============================================================================
  // Compute rate.
  //============================================================================

  return
    d_prefactor * ( d_Lnu / 1.e51 ) *
      (
         d_Enu_MeV - 2. * d_Delta + 1.2 * gsl_pow_2( d_Delta ) / d_Enu_MeV
      ) /
      gsl_pow_2( d_R6 );

}
    
//##############################################################################
// nu_n_capture_function()
//##############################################################################

double
nu_n_capture_function(
  Libnucnet__Reaction *p_reaction,
  double d_t9,
  nnt::Zone& zone
)
{

  double d_Lnu, d_Enu_MeV, d_prefactor, d_Delta, d_R6; 
  const char *s_tmp;

  if( d_t9 <= 0 )
  {
    fprintf( stderr, "Must have t9 > 0 for this function.\n" );
    exit( EXIT_FAILURE );
  }

  //============================================================================
  // Neutrino parameters.
  //============================================================================

  d_Lnu = boost::lexical_cast<double>( zone.getProperty( S_L_NU_E ) );

  d_Enu_MeV =
    3.15 * boost::lexical_cast<double>( zone.getProperty( S_T_NU_E ) );

  d_R6 =
    boost::lexical_cast<double>( zone.getProperty( nnt::s_RADIUS ) ) / 1.e6;

  //============================================================================
  // Rate parameters.
  //============================================================================

  if(
    ( s_tmp =
        Libnucnet__Reaction__getUserRateFunctionProperty(
          p_reaction,
          S_PREFACTOR,
          NULL,
          NULL
        )
    )
  )
    d_prefactor = atof( s_tmp );
  else
  {
    fprintf( stderr, " prefactor not provided.\n" );
    exit( EXIT_FAILURE );
  }

  if(
    ( s_tmp =
        Libnucnet__Reaction__getUserRateFunctionProperty(
          p_reaction,
          S_DELTA,
          NULL,
          NULL
        )
    )
  )
    d_Delta = atof( s_tmp );
  else
  {
    fprintf( stderr, " Delta not provided.\n" );
    exit( EXIT_FAILURE );
  }

  //============================================================================
  // Compute rate.
  //============================================================================

  return
    d_prefactor * ( d_Lnu / 1.e51 ) *
      (
         d_Enu_MeV + 2. * d_Delta + 1.2 * gsl_pow_2( d_Delta ) / d_Enu_MeV
      ) /
      gsl_pow_2( d_R6 );

}
    
//##############################################################################
// get_reactant_neutrino()
//##############################################################################

int
get_reactant_neutrino(
  Libnucnet__Reaction__Element *p_reactant,
  char *s_neutrino
)
{

  if( Libnucnet__Reaction__Element__isNuclide( p_reactant ) ) return 1;

  strcpy( s_neutrino, Libnucnet__Reaction__Element__getName( p_reactant ) );

  return 0;

} 

//##############################################################################
// compute_neutrino_rate()
//##############################################################################

double
compute_neutrino_rate(
  nnt::Zone& zone,
  Libnucnet__Reaction * p_reaction,
  const char *s_tnu,
  const char *s_lnu
)
{

  double d_Lnu_MeV, d_Enu_MeV, d_rate;

  d_Lnu_MeV =
    boost::lexical_cast<double>( zone.getProperty( s_lnu ) ) /
    ( GSL_CONST_CGSM_ELECTRON_VOLT * GSL_CONST_NUM_MEGA );

  d_Enu_MeV =
    3.15 * boost::lexical_cast<double>( zone.getProperty( s_tnu  ) );

  d_rate =
    compute_neutrino_cross_section(
      p_reaction,
      boost::lexical_cast<double>( zone.getProperty( s_tnu  ) )
    ) *
    ( d_Lnu_MeV / d_Enu_MeV ) /
    (
      4. *
      M_PI *
      gsl_pow_2(
        boost::lexical_cast<double>( zone.getProperty( nnt::s_RADIUS ) )
      )
    );

  return d_rate;
    
}

//##############################################################################
// compute_neutrino_cross_section().
//##############################################################################

double
compute_neutrino_cross_section(
  Libnucnet__Reaction * p_reaction,
  double d_tnu
)
{

  boost::unordered_map<std::string,NeutrinoQuantity>::iterator
    it = my_xsecs.find( Libnucnet__Reaction__getString( p_reaction ) );

  if( it != my_xsecs.end() )
  {
    return pow( 10., it->second.computeValue( d_tnu ) );
  }
  else
  {
    NeutrinoQuantity my_log10_xsec = NeutrinoQuantity( p_reaction );
    return pow( 10., my_log10_xsec.computeValue( d_tnu ) );
  }

}

//##############################################################################
// swap_neutrinos ().  Here we swap anti-neutrino_e and anti_neutrino_tau.
//##############################################################################

void
swap_neutrinos( nnt::Zone& zone )
{

  static int i_flag = 0;
  std::string s_swap;

  if( 
    (
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) )
      <
      boost::lexical_cast<double>( zone.getProperty( S_RHO_RES ) )
    ) && i_flag == 0
  )
  {

    s_swap = zone.getProperty( S_T_NU_E_BAR );

    zone.updateProperty(
      S_T_NU_E_BAR,
      zone.getProperty( S_T_NU_TAU_BAR )
    );

    zone.updateProperty(
      S_T_NU_TAU_BAR,
      s_swap
    );

    s_swap = zone.getProperty( S_L_NU_E_BAR );

    zone.updateProperty(
      S_L_NU_E,
      zone.getProperty( S_L_NU_TAU_BAR )
    );

    zone.updateProperty(
      S_L_NU_TAU_BAR,
      s_swap
    );

    i_flag = 1;

  }

}

} // namespace user
