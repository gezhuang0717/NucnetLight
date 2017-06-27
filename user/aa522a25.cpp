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
//////////////////////////////////////////////////////////////////////////////*/

////////////////////////////////////////////////////////////////////////////////
//!
//! \file
//! \brief Code for computation of approximate weak rates according to
//!        Arcones et al., A&A 522, A25 (2010) 
//!
////////////////////////////////////////////////////////////////////////////////

#include "aa522a25.h"

using namespace std;

/**
 * @brief The NucNet Tools namespace.
 */
namespace user
{

//##############################################################################
// aa522a25__update_net().
//##############################################################################

/**
 * \brief Set all the valid weak reactions in a network to use the Arcones
 *        et al. approximate weak rates.
 * \param p_my_nucnet A pointer to a Libnucnet__Net structure.
 * \param s_nuc_xpath An XPath expression giving the network subset to be
 *        updated (optional).
 */

void
aa522a25__update_net(
  Libnucnet__Net *p_net,
  const char * s_nuc_xpath
)
{

  Libnucnet__NucView * p_nuc_view =
    Libnucnet__NucView__new(
      Libnucnet__Net__getNuc( p_net ),
      s_nuc_xpath
    );

  Libnucnet__Nuc__iterateSpecies(
    Libnucnet__NucView__getNuc( p_nuc_view ),
    (Libnucnet__Species__iterateFunction)
       aa522a25__update_reactions,
    p_net
  );

  Libnucnet__NucView__free( p_nuc_view );

}

void
aa522a25__update_net(
  Libnucnet__Net *p_net
)
{

  Libnucnet__Nuc__iterateSpecies(
    Libnucnet__Net__getNuc( p_net ),
    (Libnucnet__Species__iterateFunction)
       aa522a25__update_reactions,
    p_net
  );

}

//##############################################################################
// aa522a25__update_reactions().
//##############################################################################

int
aa522a25__update_reactions(
  Libnucnet__Species *p_species,
  Libnucnet__Net *p_net
)
{

  Libnucnet__Species *p_product;
  Libnucnet__Reaction *p_reaction;
  char s_b[32];

  //============================================================================
  // Set aa522a25 B factor.
  //============================================================================

  if( Libnucnet__Species__getA( p_species ) == 1 )
    strcpy( s_b, "5.76" );
  else 
    strcpy( s_b, "4.6" );

  //============================================================================
  // Check that electron capture or beta+ can happen.
  //============================================================================

  if(
    Libnucnet__Species__getZ( p_species ) > 1
    ||
    (
      Libnucnet__Species__getZ( p_species ) == 1 && 
      Libnucnet__Species__getA( p_species ) < 2
    )
  )
  {

    //--------------------------------------------------------------------------
    // Electron capture.
    //--------------------------------------------------------------------------

    p_product =
      Libnucnet__Nuc__getSpeciesByZA(
        Libnucnet__Net__getNuc( p_net ),
        Libnucnet__Species__getZ( p_species ) - 1,
        Libnucnet__Species__getA( p_species ),
        NULL
      );

    if( p_product )
    {

      p_reaction = Libnucnet__Reaction__new(); 

      Libnucnet__Reaction__updateSource(
        p_reaction,
        "aa522a25"
      );

      Libnucnet__Reaction__addReactant(
        p_reaction,
        Libnucnet__Species__getName( p_species )
      );

      Libnucnet__Reaction__addReactant( p_reaction, nnt::s_ELECTRON );

      Libnucnet__Reaction__addProduct(
        p_reaction,
        Libnucnet__Species__getName( p_product )
      );

      Libnucnet__Reaction__addProduct( p_reaction, nnt::s_NEUTRINO_E );

      Libnucnet__Reaction__setUserRateFunctionKey(
        p_reaction,
        s_AA522A25_ELECTRON_CAPTURE
      );

      Libnucnet__Reac__removeReaction(
        Libnucnet__Net__getReac( p_net ),
        p_reaction
      );

      Libnucnet__Reaction__updateUserRateFunctionProperty(
        p_reaction,
        s_AA522A25_B_FACTOR,
        NULL,
        NULL,
        s_b
      );

      Libnucnet__Reac__addReaction(
        Libnucnet__Net__getReac( p_net ),
        p_reaction
      );

    }

    //--------------------------------------------------------------------------
    // Beta plus.
    //--------------------------------------------------------------------------

    p_product =
      Libnucnet__Nuc__getSpeciesByZA(
        Libnucnet__Net__getNuc( p_net ),
        Libnucnet__Species__getZ( p_species ) - 1,
        Libnucnet__Species__getA( p_species ),
        NULL
      );

    if( p_product )
    {
 
      p_reaction = Libnucnet__Reaction__new(); 

      Libnucnet__Reaction__updateSource(
        p_reaction,
        "aa522a25"
      );

      Libnucnet__Reaction__addReactant(
        p_reaction,
        Libnucnet__Species__getName( p_species )
      );

      Libnucnet__Reaction__addProduct(
        p_reaction,
        Libnucnet__Species__getName( p_product )
      );

      Libnucnet__Reaction__addProduct( p_reaction, nnt::s_POSITRON );

      Libnucnet__Reaction__addProduct( p_reaction, nnt::s_NEUTRINO_E );

      Libnucnet__Reaction__setUserRateFunctionKey(
        p_reaction,
        s_AA522A25_BETA_PLUS
      );

      Libnucnet__Reac__removeReaction(
        Libnucnet__Net__getReac( p_net ),
        p_reaction
      );

      Libnucnet__Reaction__updateUserRateFunctionProperty(
        p_reaction,
        s_AA522A25_B_FACTOR,
        NULL,
        NULL,
        s_b
      );

      Libnucnet__Reac__addReaction(
        Libnucnet__Net__getReac( p_net ),
        p_reaction
      );

    }

  }

  //============================================================================
  // Beta - decay.
  //============================================================================

  p_product =
    Libnucnet__Nuc__getSpeciesByZA(
      Libnucnet__Net__getNuc( p_net ),
      Libnucnet__Species__getZ( p_species ) + 1,
      Libnucnet__Species__getA( p_species ),
      NULL
    );

  if( p_product )
  {

    p_reaction = Libnucnet__Reaction__new(); 

    Libnucnet__Reaction__updateSource(
      p_reaction,
      "aa522a25"
    );

    Libnucnet__Reaction__addReactant(
      p_reaction,
      Libnucnet__Species__getName( p_species )
    );

    Libnucnet__Reaction__addProduct(
      p_reaction,
      Libnucnet__Species__getName( p_product )
    );

    Libnucnet__Reaction__addProduct( p_reaction, nnt::s_ELECTRON );

    Libnucnet__Reaction__addProduct( p_reaction, nnt::s_ANTI_NEUTRINO_E );

    Libnucnet__Reaction__setUserRateFunctionKey(
      p_reaction,
      s_AA522A25_BETA_MINUS
    );

    Libnucnet__Reac__removeReaction(
      Libnucnet__Net__getReac( p_net ),
      p_reaction
    );

    Libnucnet__Reaction__updateUserRateFunctionProperty(
      p_reaction,
      s_AA522A25_B_FACTOR,
      NULL,
      NULL,
      s_b
    );

    Libnucnet__Reac__addReaction(
      Libnucnet__Net__getReac( p_net ),
      p_reaction
    );

  }

  //============================================================================
  // Positron capture.
  //============================================================================

  p_product =
    Libnucnet__Nuc__getSpeciesByZA(
      Libnucnet__Net__getNuc( p_net ),
      Libnucnet__Species__getZ( p_species ) + 1,
      Libnucnet__Species__getA( p_species ),
      NULL
    );

  if( p_product )
  {

    p_reaction = Libnucnet__Reaction__new(); 
  
    Libnucnet__Reaction__updateSource(
      p_reaction,
      "aa522a25"
    );

    Libnucnet__Reaction__addReactant(
      p_reaction,
      Libnucnet__Species__getName( p_species )
    );

    Libnucnet__Reaction__addReactant( p_reaction, nnt::s_POSITRON );

    Libnucnet__Reaction__addProduct(
      p_reaction,
      Libnucnet__Species__getName( p_product )
    );

    Libnucnet__Reaction__addProduct( p_reaction, nnt::s_ANTI_NEUTRINO_E );

    Libnucnet__Reaction__setUserRateFunctionKey(
      p_reaction,
      s_AA522A25_POSITRON_CAPTURE
    );

    Libnucnet__Reac__removeReaction(
      Libnucnet__Net__getReac( p_net ),
      p_reaction
    );

    Libnucnet__Reaction__updateUserRateFunctionProperty(
      p_reaction,
      s_AA522A25_B_FACTOR,
      NULL,
      NULL,
      s_b
    );

    Libnucnet__Reac__addReaction(
      Libnucnet__Net__getReac( p_net ),
      p_reaction
    );

  }

  return 1;

}

//##############################################################################
// aa522a25__compute_rate().
//##############################################################################

double
aa522a25__compute_rate(
  Libnucnet__Reaction * p_reaction,
  Libnucnet__Net * p_net,
  double d_electron_mass,
  double d_t9,
  double d_eta_F,
  double d_mu_nue_kT
)
{

  double d_result;
     
  if(
    strcmp(
      Libnucnet__Reaction__getRateFunctionKey( p_reaction ),
      s_AA522A25_ELECTRON_CAPTURE
    ) == 0
  )
  {
    d_result =
      aa522a25__compute_electron_capture_integral(
        p_reaction,
        p_net,
        d_electron_mass,
        d_t9,
        d_eta_F,
        d_mu_nue_kT,
        2
      );
  }
  else if(
    strcmp(
      Libnucnet__Reaction__getRateFunctionKey( p_reaction ),
      s_AA522A25_BETA_PLUS
    ) == 0
  )
  {
    d_result =
      aa522a25__compute_beta_plus_integral(
        p_reaction,
        p_net,
        d_electron_mass,
        d_t9,
        d_eta_F,
        d_mu_nue_kT,
        2
      );
  }
  else if(
    strcmp(
      Libnucnet__Reaction__getRateFunctionKey( p_reaction ),
      s_AA522A25_POSITRON_CAPTURE
    ) == 0
  )
  {
    d_result =
      aa522a25__compute_positron_capture_integral(
        p_reaction,
        p_net,
        d_electron_mass,
        d_t9,
        d_eta_F,
        d_mu_nue_kT,
        2
      );
  }
  else if(
    strcmp(
      Libnucnet__Reaction__getRateFunctionKey( p_reaction ),
      s_AA522A25_BETA_MINUS
    ) == 0
  )
  {
    d_result =
      aa522a25__compute_beta_minus_integral(
        p_reaction,
        p_net,
        d_electron_mass,
        d_t9,
        d_eta_F,
        d_mu_nue_kT,
        2
      );
  }
  else
  {
    std::cerr <<
      "No such rate for AA522A25 reaction rate function." << std::endl;
    exit( EXIT_FAILURE );
  }

  if(
      !Libnucnet__Reaction__getUserRateFunctionProperty(
        p_reaction,
        s_AA522A25_B_FACTOR,
        NULL,
        NULL
      )
  )
  {
    std::cerr <<
      "No such user rate property." << std::endl;
    exit( EXIT_FAILURE );
  }

  return
    d_result *
    atof(
      Libnucnet__Reaction__getUserRateFunctionProperty(
        p_reaction,
        s_AA522A25_B_FACTOR,
        NULL,
        NULL
      )
    ) *
    M_LN2 /
    (
      gsl_pow_5( d_electron_mass )
      *
      d_AA522A25_WEAK_K
    );


}

//##############################################################################
// aa522a25__compute_electron_capture_integral().
//##############################################################################

double
aa522a25__compute_electron_capture_integral(
  Libnucnet__Reaction * p_reaction,
  Libnucnet__Net * p_net,
  double d_electron_mass,
  double d_t9,
  double d_eta_F,
  double d_mu_nue_kT,
  int i_neutrino_exponent
)
{

  double d_result = 0, d_error, d_lower_limit, d_upper_limit;
  gsl_integration_workspace * p_work_space;
  gsl_function F;
  aa522a25_function_data function_data;
     
  //============================================================================
  // Set up function data structure.
  //============================================================================

  function_data.dEtaF = d_eta_F;
  function_data.dMunuekT = d_mu_nue_kT;
  function_data.dkT = nnt::compute_kT_in_MeV( d_t9 ),
  function_data.dQ =
    nnt::compute_reaction_nuclear_Qvalue(
      p_net,
      p_reaction,
      d_electron_mass
    );
  function_data.iNeutrinoExponent = i_neutrino_exponent;
    
  F.params = &function_data;

  F.function = &aa522a25__electron_capture_integrand;

  function_data.dQ *= -1;

  //============================================================================
  // Integration limits.
  //============================================================================

  d_lower_limit =
    GSL_MAX_DBL(
      function_data.dQ,
      d_electron_mass
    );

  d_upper_limit = function_data.dkT * ( 700. + d_eta_F );

  if( d_upper_limit < d_lower_limit ) return 0.;

  //============================================================================
  // Set up workspace and integrate.
  //============================================================================

  p_work_space =
    gsl_integration_workspace_alloc( (size_t) i_AA522A25_WORKSPACE );

  gsl_integration_qags(
    &F,
    d_lower_limit,
    d_upper_limit,
    d_AA522A25_EPSABS,
    d_AA522A25_EPSREL,
    (size_t) i_AA522A25_LIMIT,
    p_work_space,
    &d_result,
    &d_error
  );

  gsl_integration_workspace_free( p_work_space );

  return d_result;

}

//##############################################################################
// aa522a25__compute_beta_plus_integral().
//##############################################################################

double
aa522a25__compute_beta_plus_integral(
  Libnucnet__Reaction * p_reaction,
  Libnucnet__Net * p_net,
  double d_electron_mass,
  double d_t9,
  double d_eta_F,
  double d_mu_nue_kT,
  int i_neutrino_exponent
)
{

  double d_result = 0, d_error;
  gsl_integration_workspace * p_work_space;
  gsl_function F;
  aa522a25_function_data function_data;
     
  //============================================================================
  // Set up workspace and function data structure.
  //============================================================================

  p_work_space =
    gsl_integration_workspace_alloc( (size_t) i_AA522A25_WORKSPACE );

  function_data.dEtaF = d_eta_F;
  function_data.dMunuekT = d_mu_nue_kT;
  function_data.dkT = nnt::compute_kT_in_MeV( d_t9 ),
  function_data.dQ =
    nnt::compute_reaction_nuclear_Qvalue(
      p_net,
      p_reaction,
      d_electron_mass
    );
  function_data.iNeutrinoExponent = i_neutrino_exponent;
    
  F.function = &aa522a25__beta_plus_integrand;
  F.params = &function_data;

  function_data.dQ *= -1;

  //============================================================================
  // Integration.
  //============================================================================

  if( 
    d_electron_mass >
    -function_data.dQ
  )
    d_result = 0.;
  else
    gsl_integration_qags(
      &F,
      d_electron_mass,
      -function_data.dQ,
      d_AA522A25_EPSABS,
      d_AA522A25_EPSREL,
      (size_t) i_AA522A25_LIMIT,
      p_work_space,
      &d_result,
      &d_error
    );

  gsl_integration_workspace_free( p_work_space );

  return d_result;

}

//##############################################################################
// aa522a25__compute_position_capture_integral().
//##############################################################################

double
aa522a25__compute_positron_capture_integral(
  Libnucnet__Reaction * p_reaction,
  Libnucnet__Net * p_net,
  double d_electron_mass,
  double d_t9,
  double d_eta_F,
  double d_mu_nue_kT,
  int i_neutrino_exponent
)
{

  double d_result = 0, d_error;
  gsl_integration_workspace * p_work_space;
  gsl_function F;
  aa522a25_function_data function_data;
     
  //============================================================================
  // Set up workspace and function data structure.
  //============================================================================

  p_work_space =
    gsl_integration_workspace_alloc( (size_t) i_AA522A25_WORKSPACE );

  function_data.dEtaF = d_eta_F;
  function_data.dMunuekT = d_mu_nue_kT;
  function_data.dkT = nnt::compute_kT_in_MeV( d_t9 ),
  function_data.dQ =
    nnt::compute_reaction_nuclear_Qvalue(
      p_net,
      p_reaction,
      d_electron_mass
    );
  function_data.iNeutrinoExponent = i_neutrino_exponent;
    
  F.function = &aa522a25__positron_capture_integrand;
  F.params = &function_data;

  //============================================================================
  // Integration.
  //============================================================================

  gsl_integration_qagiu(
    &F,
    GSL_MAX_DBL(
      -function_data.dQ,
      d_electron_mass
    ),
    d_AA522A25_EPSABS,
    d_AA522A25_EPSREL,
    (size_t) i_AA522A25_LIMIT,
    p_work_space,
    &d_result,
    &d_error
  );

  gsl_integration_workspace_free( p_work_space );

  return d_result;

}

//##############################################################################
// aa522a25__compute_beta_minus_integral().
//##############################################################################

double
aa522a25__compute_beta_minus_integral(
  Libnucnet__Reaction * p_reaction,
  Libnucnet__Net * p_net,
  double d_electron_mass,
  double d_t9,
  double d_eta_F,
  double d_mu_nue_kT,
  int i_neutrino_exponent
)
{

  double d_result = 0, d_error;
  gsl_integration_workspace * p_work_space;
  gsl_function F;
  aa522a25_function_data function_data;
     
  //============================================================================
  // Set up workspace and function data structure.
  //============================================================================

  p_work_space =
    gsl_integration_workspace_alloc( (size_t) i_AA522A25_WORKSPACE );

  function_data.dEtaF = d_eta_F;
  function_data.dMunuekT = d_mu_nue_kT;
  function_data.dkT = nnt::compute_kT_in_MeV( d_t9 ),
  function_data.dQ =
    nnt::compute_reaction_nuclear_Qvalue(
      p_net,
      p_reaction,
      d_electron_mass
    );
  function_data.iNeutrinoExponent = i_neutrino_exponent;
    
  F.function = &aa522a25__beta_minus_integrand;
  F.params = &function_data;

  if( 
    d_electron_mass > function_data.dQ
  )
    d_result = 0.;
  else
    gsl_integration_qags(
      &F,
      d_electron_mass,
      function_data.dQ,
      d_AA522A25_EPSABS,
      d_AA522A25_EPSREL,
      (size_t) i_AA522A25_LIMIT,
      p_work_space,
      &d_result,
      &d_error
    );

  gsl_integration_workspace_free( p_work_space );

  return d_result;

}

//##############################################################################
// aa522a25__electron_capture_integrand().
//##############################################################################

double
aa522a25__electron_capture_integrand(
  double d_E,
  void *p_data
)
{

  double d_result;

  aa522a25_function_data *p_function_data = ( aa522a25_function_data * ) p_data;

  d_result =
    pow( d_E - p_function_data->dQ, p_function_data->iNeutrinoExponent );

  d_result *=
    gsl_pow_2( d_E ) *
    nnt::compute_fermi_dirac_factor(
      d_E / p_function_data->dkT,
      p_function_data->dEtaF
    );

  if( p_function_data->dMunuekT == GSL_NEGINF )
    return d_result;
  else
    return
      d_result *
      nnt::compute_one_minus_fermi_dirac_factor(
        d_E / p_function_data->dkT - p_function_data->dQ / p_function_data->dkT,
        p_function_data->dMunuekT
      );

}

//##############################################################################
// aa522a25__beta_plus_integrand().
//##############################################################################

double
aa522a25__beta_plus_integrand(
  double d_E,
  void *p_data
)
{

  double d_result;

  aa522a25_function_data *p_function_data = ( aa522a25_function_data * ) p_data;

  d_result =
    pow( -p_function_data->dQ - d_E, p_function_data->iNeutrinoExponent );

  d_result *=
    gsl_pow_2( d_E ) *
    nnt::compute_one_minus_fermi_dirac_factor(
      d_E / p_function_data->dkT,
      -p_function_data->dEtaF
    );

  if( p_function_data->dMunuekT == GSL_NEGINF )
    return d_result;
  else
    return
      d_result *
      nnt::compute_one_minus_fermi_dirac_factor(
        -p_function_data->dQ / p_function_data->dkT -
           d_E / p_function_data->dkT,
        p_function_data->dMunuekT
      );

}

//##############################################################################
// aa522a25__positron_capture_integrand().
//##############################################################################

double
aa522a25__positron_capture_integrand(
  double d_E,
  void *p_data
)
{

  double d_result;

  aa522a25_function_data *p_function_data = ( aa522a25_function_data * ) p_data;

  d_result =
    pow( d_E + p_function_data->dQ, p_function_data->iNeutrinoExponent );

  d_result *=
    gsl_pow_2( d_E ) *
    nnt::compute_fermi_dirac_factor(
      d_E / p_function_data->dkT,
      -p_function_data->dEtaF
    );

  if( p_function_data->dMunuekT == GSL_NEGINF )
    return d_result;
  else
    return
      d_result *
      nnt::compute_one_minus_fermi_dirac_factor(
        d_E / p_function_data->dkT + p_function_data->dQ / p_function_data->dkT,
        -p_function_data->dMunuekT
      );

}

//##############################################################################
// aa522a25__beta_minus_integrand().
//##############################################################################

double
aa522a25__beta_minus_integrand(
  double d_E,
  void *p_data
)
{

  double d_result;

  aa522a25_function_data *p_function_data = ( aa522a25_function_data * ) p_data;

  d_result =
    pow( p_function_data->dQ - d_E, p_function_data->iNeutrinoExponent );

  d_result *=
    gsl_pow_2( d_E ) *
    nnt::compute_one_minus_fermi_dirac_factor(
      d_E / p_function_data->dkT,
      p_function_data->dEtaF
    );

  if( p_function_data->dMunuekT == GSL_NEGINF )
    return d_result;
  else
    return
      d_result *
      nnt::compute_one_minus_fermi_dirac_factor(
        p_function_data->dQ / p_function_data->dkT - d_E / p_function_data->dkT,
        -p_function_data->dMunuekT
      );

}

//##############################################################################
// aa522a25__compute_reaction_neutrino_energy_loss_rate().
//##############################################################################

double
aa522a25__compute_reaction_neutrino_energy_loss_rate(
  Libnucnet__Reaction * p_reaction,
  Libnucnet__Net * p_net,
  double d_electron_mass,
  double d_t9,
  double d_eta_F,
  double d_mu_nue_kT
)
{

  double d_result = 0, d_error, d_limit;
  gsl_integration_workspace * p_work_space =
    gsl_integration_workspace_alloc( (size_t) i_AA522A25_WORKSPACE );
  gsl_function F;
  aa522a25_function_data function_data;

  function_data.dkT = nnt::compute_kT_in_MeV( d_t9 );

  function_data.dEtaF = d_eta_F;

  function_data.dMunuekT = d_mu_nue_kT;

  function_data.iNeutrinoExponent = 3;

  function_data.dQ =
    nnt::compute_reaction_nuclear_Qvalue(
      p_net,
      p_reaction,
      d_electron_mass
    );

  F.params = &function_data;

  if(
    strcmp(
      Libnucnet__Reaction__getRateFunctionKey( p_reaction ),
      s_AA522A25_ELECTRON_CAPTURE
    ) == 0
  )
  {

    F.function = &aa522a25__electron_capture_integrand;

    function_data.dQ *= -1;

    d_limit =
      GSL_MAX_DBL(
        function_data.dQ,
        d_electron_mass
      );

    gsl_integration_qagiu(
      &F,
      d_limit,
      d_AA522A25_EPSABS,
      d_AA522A25_EPSREL,
      (size_t) i_AA522A25_LIMIT,
      p_work_space,
      &d_result,
      &d_error
    );

  }
  else if(
    strcmp(
      Libnucnet__Reaction__getRateFunctionKey( p_reaction ),
      s_AA522A25_BETA_PLUS
    ) == 0
  )
  {

    F.function = &aa522a25__beta_plus_integrand;

    function_data.dQ *= -1;

    if( 
      d_electron_mass > -function_data.dQ
    )
      d_result = 0.;
    else
      gsl_integration_qags(
        &F,
        d_electron_mass,
        -function_data.dQ,
        d_AA522A25_EPSABS,
        d_AA522A25_EPSREL,
        (size_t) i_AA522A25_LIMIT,
        p_work_space,
        &d_result,
        &d_error
      );

  }
  else if(
    strcmp(
      Libnucnet__Reaction__getRateFunctionKey( p_reaction ),
      s_AA522A25_POSITRON_CAPTURE
    ) == 0
  )
  {

    F.function = &aa522a25__positron_capture_integrand;

    gsl_integration_qagiu(
      &F,
      GSL_MAX_DBL(
        -function_data.dQ,
        d_electron_mass
      ),
      d_AA522A25_EPSABS,
      d_AA522A25_EPSREL,
      (size_t) i_AA522A25_LIMIT,
      p_work_space,
      &d_result,
      &d_error
    );

  }
  else if(
    strcmp(
      Libnucnet__Reaction__getRateFunctionKey( p_reaction ),
      s_AA522A25_BETA_MINUS
    ) == 0
  )
  {

    F.function = &aa522a25__beta_minus_integrand;

    if( 
      d_electron_mass >
      function_data.dQ
    )
      d_result = 0.;
    else
      gsl_integration_qags(
        &F,
        d_electron_mass,
        function_data.dQ,
        d_AA522A25_EPSABS,
        d_AA522A25_EPSREL,
        (size_t) i_AA522A25_LIMIT,
        p_work_space,
        &d_result,
        &d_error
      );

  }
  else
  {
    std::cerr <<
      "No such rate for AA522A25 reaction rate function." << std::endl;
    exit( EXIT_FAILURE );
  }

  gsl_integration_workspace_free( p_work_space );

  if(
      !Libnucnet__Reaction__getUserRateFunctionProperty(
        p_reaction,
        s_AA522A25_B_FACTOR,
        NULL,
        NULL
      )
  )
  {
    std::cerr <<
      "No such user rate property." << std::endl;
    exit( EXIT_FAILURE );
  }

  return
    d_result *
    atof(
      Libnucnet__Reaction__getUserRateFunctionProperty(
        p_reaction,
        s_AA522A25_B_FACTOR,
        NULL,
        NULL
      )
    ) *
    M_LN2 /
    (
      gsl_pow_5( d_electron_mass )
      *
      d_AA522A25_WEAK_K
    );


}

//##############################################################################
// aa522a25__compute_reaction_average_neutrino_energy().
//##############################################################################

double
aa522a25__compute_reaction_average_neutrino_energy(
  Libnucnet__Reaction * p_reaction,
  Libnucnet__Net * p_net,
  double d_electron_mass,
  double d_t9,
  double d_eta_F,
  double d_mu_nue_kT
)
{

  double d_denom;

  if( nnt::is_electron_capture_reaction( p_reaction ) )
  {

    d_denom = 
      aa522a25__compute_electron_capture_integral(
        p_reaction,
        p_net,
        d_electron_mass,
        d_t9,
        d_eta_F,
        d_mu_nue_kT,
        2
      );

    if( GSL_SIGN( d_denom ) == GSL_SIGN( -d_denom ) )
      return 0;
    else
      return
        aa522a25__compute_electron_capture_integral(
          p_reaction,
          p_net,
          d_electron_mass,
          d_t9,
          d_eta_F,
          d_mu_nue_kT,
          3
        ) / d_denom;

  }
  else if( nnt::is_beta_minus_reaction( p_reaction ) )
  {

    d_denom = 
      aa522a25__compute_beta_minus_integral(
        p_reaction,
        p_net,
        d_electron_mass,
        d_t9,
        d_eta_F,
        d_mu_nue_kT,
        2
      );

    if( GSL_SIGN( d_denom ) == GSL_SIGN( -d_denom ) )
      return 0;
    else
      return
        aa522a25__compute_beta_minus_integral(
          p_reaction,
          p_net,
          d_electron_mass,
          d_t9,
          d_eta_F,
          d_mu_nue_kT,
          3
        ) / d_denom;

  }
  else if( nnt::is_beta_plus_reaction( p_reaction ) )
  {

    d_denom = 
      aa522a25__compute_beta_plus_integral(
        p_reaction,
        p_net,
        d_electron_mass,
        d_t9,
        d_eta_F,
        d_mu_nue_kT,
        2
      );

    if( GSL_SIGN( d_denom ) == GSL_SIGN( -d_denom ) )
      return 0;
    else
      return
        aa522a25__compute_beta_plus_integral(
          p_reaction,
          p_net,
          d_electron_mass,
          d_t9,
          d_eta_F,
          d_mu_nue_kT,
          3
        ) / d_denom;

  }
  else if( nnt::is_positron_capture_reaction( p_reaction ) )
  {

    d_denom = 
      aa522a25__compute_positron_capture_integral(
        p_reaction,
        p_net,
        d_electron_mass,
        d_t9,
        d_eta_F,
        d_mu_nue_kT,
        2
      );

    if( GSL_SIGN( d_denom ) == GSL_SIGN( -d_denom ) )
      return 0;
    else
      return
        aa522a25__compute_positron_capture_integral(
          p_reaction,
          p_net,
          d_electron_mass,
          d_t9,
          d_eta_F,
          d_mu_nue_kT,
          3
        ) / d_denom;

  }
  else
    return 0;

}
     
} // namespace user
