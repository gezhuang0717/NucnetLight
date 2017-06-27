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
//! \brief Code for user-defined rate functions.
//!
////////////////////////////////////////////////////////////////////////////////

#include "user_rate_functions.h"

/**
 * @brief A namespace for user-defined rate functions.
 */
namespace user
{

//##############################################################################
// register_my_rate_functions().
//##############################################################################

void
register_my_rate_functions( Libnucnet__Reac *p_reac )
{

  //============================================================================
  // Register two-d weak rates and set deallocator.
  //============================================================================

  Libnucnet__Reac__registerUserRateFunction(
    p_reac,
    nnt::s_TWO_D_WEAK_RATES, 
    (Libnucnet__Reaction__userRateFunction)
       compute_two_d_weak_rate
  );

  Libnucnet__Reac__setUserRateFunctionDataDeallocator(
    p_reac,
    nnt::s_TWO_D_WEAK_RATES,
    (Libnucnet__Reaction__user_rate_function_data_deallocator) free
  );

  //============================================================================
  // Register two-d weak rates log10 ft and set the deallocator.
  //============================================================================

  Libnucnet__Reac__registerUserRateFunction(
    p_reac,
    nnt::s_TWO_D_WEAK_RATES_LOG10_FT, 
    (Libnucnet__Reaction__userRateFunction)
       my_log10_ft_rate_function
  );

  Libnucnet__Reac__setUserRateFunctionDataDeallocator(
    p_reac,
    nnt::s_TWO_D_WEAK_RATES_LOG10_FT, 
    (Libnucnet__Reaction__user_rate_function_data_deallocator)
       free
  );

  //============================================================================
  // Register approximate weak rate functions (A&A vol. 522, A25) and set their
  // deallocators.
  //============================================================================

  Libnucnet__Reac__registerUserRateFunction(
    p_reac,
    s_AA522A25_ELECTRON_CAPTURE,
    (Libnucnet__Reaction__userRateFunction)
       my_approximate_weak_rate_function
  );

  Libnucnet__Reac__setUserRateFunctionDataDeallocator(
    p_reac,
    s_AA522A25_ELECTRON_CAPTURE,
    (Libnucnet__Reaction__user_rate_function_data_deallocator) free
  );

  Libnucnet__Reac__registerUserRateFunction(
    p_reac,
    s_AA522A25_BETA_PLUS,
    (Libnucnet__Reaction__userRateFunction)
       my_approximate_weak_rate_function
  );

  Libnucnet__Reac__setUserRateFunctionDataDeallocator(
    p_reac,
    s_AA522A25_BETA_PLUS,
    (Libnucnet__Reaction__user_rate_function_data_deallocator) free
  );

  Libnucnet__Reac__registerUserRateFunction(
    p_reac,
    s_AA522A25_POSITRON_CAPTURE,
    (Libnucnet__Reaction__userRateFunction)
       my_approximate_weak_rate_function
  );

  Libnucnet__Reac__setUserRateFunctionDataDeallocator(
    p_reac,
    s_AA522A25_POSITRON_CAPTURE,
    (Libnucnet__Reaction__user_rate_function_data_deallocator) free
  );

  Libnucnet__Reac__registerUserRateFunction(
    p_reac,
    s_AA522A25_BETA_MINUS,
    (Libnucnet__Reaction__userRateFunction)
       my_approximate_weak_rate_function
  );

  Libnucnet__Reac__setUserRateFunctionDataDeallocator(
    p_reac,
    s_AA522A25_BETA_MINUS,
    (Libnucnet__Reaction__user_rate_function_data_deallocator) free
  );

  //============================================================================
  // Register neutrino rate functions.
  //============================================================================

  register_my_neutrino_rate_functions( p_reac );

}

//##############################################################################
// update_my_two_d_weak_rate_functions_data().
//##############################################################################

void
update_my_two_d_weak_rate_functions_data(
  Libnucnet__Zone *p_zone,
  double d_electron_mass,
  double d_rhoe,
  double d_eta_F,
  double d_mu_nue_kT
)
{

  typedef struct {
    Libnucnet__Net *pNet;
    double dElectronMass;
    double dRhoe;
    double dEtaF;
    double dMuNuekT;
  } work;

  work *p_work;
  double *p_rhoe;

  //============================================================================
  // Set data for two-d weak rates. 
  //============================================================================

  p_rhoe = ( double * ) malloc( sizeof( double ) );

  *p_rhoe = d_rhoe;

  Libnucnet__Zone__updateDataForUserRateFunction(
    p_zone,
    nnt::s_TWO_D_WEAK_RATES,
    p_rhoe
  );

  //============================================================================
  // Set data for two-d weak rates log10 ft. 
  //============================================================================

  p_work = ( work * ) malloc( sizeof( work ) );

  if( !p_work )
  {
    fprintf( stderr, "Couldn't allocate memory for work structure.\n" );
    exit( EXIT_FAILURE );
  }

  p_work->pNet = Libnucnet__Zone__getNet( p_zone );
  p_work->dRhoe = d_rhoe;

  p_work->dElectronMass = d_electron_mass;

  p_work->dEtaF = d_eta_F;

  p_work->dMuNuekT = d_mu_nue_kT;

  Libnucnet__Zone__updateDataForUserRateFunction(
    p_zone,
    nnt::s_TWO_D_WEAK_RATES_LOG10_FT,
    p_work
  );

}

//##############################################################################
// update_my_approximate_weak_rate_functions_data().
//##############################################################################

void
update_my_approximate_weak_rate_functions_data(
  Libnucnet__Zone *p_zone,
  double d_electron_mass,
  double d_eta_F,
  double d_mu_nue_kT

)
{

  typedef struct {
    Libnucnet__Net *pNet;
    double dElectronMassMeV;
    double dEtaF;
    double dMuNuekT;
  } work;

  work * p_work_ec, * p_work_bp, * p_work_pc, * p_work_bm;

  //============================================================================
  // Set data for approximate weak rates.
  //============================================================================
  
  p_work_ec = ( work * ) malloc( sizeof( work ) );

  if( !p_work_ec )
  {
    std::cerr << "Couldn't allocate work structure." << std::endl;
    exit( EXIT_FAILURE );
  } 

  p_work_ec->pNet = Libnucnet__Zone__getNet( p_zone );
  p_work_ec->dElectronMassMeV = d_electron_mass;
  p_work_ec->dEtaF = d_eta_F;
  p_work_ec->dMuNuekT = d_mu_nue_kT;

  //============================================================================
  // Set data for approximate electron capture reaction.
  //============================================================================

  Libnucnet__Zone__updateDataForUserRateFunction(
    p_zone,
    s_AA522A25_ELECTRON_CAPTURE,
    p_work_ec
  );

  //============================================================================
  // Set data for approximate weak rates.
  //============================================================================

  p_work_bp = ( work * ) malloc( sizeof( work ) );

  if( !p_work_bp )
  {
    std::cerr << "Couldn't allocate work structure." << std::endl;
    exit( EXIT_FAILURE );
  } 

  p_work_bp->pNet = Libnucnet__Zone__getNet( p_zone );
  p_work_bp->dElectronMassMeV = d_electron_mass;
  p_work_bp->dEtaF = d_eta_F;
  p_work_bp->dMuNuekT = d_mu_nue_kT;

  //============================================================================
  // Set data for approximate bp reaction.
  //============================================================================

  Libnucnet__Zone__updateDataForUserRateFunction(
    p_zone,
    s_AA522A25_BETA_PLUS,
    p_work_bp
  );

  //============================================================================
  // Set data for approximate weak rates.
  //============================================================================

  p_work_pc = ( work * ) malloc( sizeof( work ) );

  if( !p_work_pc )
  {
    std::cerr << "Couldn't allocate work structure." << std::endl;
    exit( EXIT_FAILURE );
  } 

  p_work_pc->pNet = Libnucnet__Zone__getNet( p_zone );
  p_work_pc->dElectronMassMeV = d_electron_mass;
  p_work_pc->dEtaF = d_eta_F;
  p_work_pc->dMuNuekT = d_mu_nue_kT;

  //============================================================================
  // Set data for approximate positron capture reaction.
  //============================================================================

  Libnucnet__Zone__updateDataForUserRateFunction(
    p_zone,
    s_AA522A25_POSITRON_CAPTURE,
    p_work_pc
  );

  //============================================================================
  // Set data for approximate weak rates.
  //============================================================================

  p_work_bm = ( work * ) malloc( sizeof( work ) );

  if( !p_work_bm )
  {
    std::cerr << "Couldn't allocate work structure." << std::endl;
    exit( EXIT_FAILURE );
  } 

  p_work_bm->pNet = Libnucnet__Zone__getNet( p_zone );
  p_work_bm->dElectronMassMeV = d_electron_mass;
  p_work_bm->dEtaF = d_eta_F;
  p_work_bm->dMuNuekT = d_mu_nue_kT;

  //============================================================================
  // Set data for approximate beta-minus reaction.
  //============================================================================

  Libnucnet__Zone__updateDataForUserRateFunction(
    p_zone,
    s_AA522A25_BETA_MINUS,
    p_work_bm
  );

}

//##############################################################################
// update_my_rate_functions_data().
//##############################################################################

void
update_my_rate_functions_data(
  nnt::Zone &zone
)
{

  Libstatmech__Fermion * p_electron;
  double d_mue_kT, d_eta_F;

  if(
    Libnucnet__Reac__isRegisteredRateFunction(
      Libnucnet__Net__getReac(
        Libnucnet__Zone__getNet( zone.getNucnetZone() )
      ),
       nnt::s_TWO_D_WEAK_RATES
    ) ||
    Libnucnet__Reac__isRegisteredRateFunction(
      Libnucnet__Net__getReac(
        Libnucnet__Zone__getNet( zone.getNucnetZone() )
      ),
      nnt::s_TWO_D_WEAK_RATES_LOG10_FT
    ) ||
    ( zone.hasProperty( nnt::s_USE_APPROXIMATE_WEAK_RATES ) &&
      zone.getProperty( nnt::s_USE_APPROXIMATE_WEAK_RATES ) == "yes"
    )
  )
  {

    p_electron = 
      Libstatmech__Fermion__new(
        nnt::s_ELECTRON,
        nnt::d_ELECTRON_MASS_IN_MEV,
        2,
        -1
      );

    d_mue_kT =
      Libstatmech__Fermion__computeChemicalPotential(
        p_electron,
        boost::lexical_cast<double>(
          zone.getProperty( nnt::s_T9 )
        ) * GSL_CONST_NUM_GIGA,
        boost::lexical_cast<double>(
          zone.getProperty( nnt::s_RHO )
        ) * Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 ) *
          GSL_CONST_NUM_AVOGADRO,
        NULL,
        NULL
      );

    d_eta_F =
      d_mue_kT
      +
      Libstatmech__Fermion__getRestMass( p_electron ) /
      nnt::compute_kT_in_MeV( 
        boost::lexical_cast<double>(
        zone.getProperty( nnt::s_T9 )
      )
    );

    update_my_two_d_weak_rate_functions_data(
      zone.getNucnetZone(),
      Libstatmech__Fermion__getRestMass( p_electron ),
      boost::lexical_cast<double>(
        zone.getProperty( nnt::s_RHO )
      ) * Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 ),
      d_eta_F,
      zone.getElectronNeutrinoChemicalPotential()
    );

    update_my_approximate_weak_rate_functions_data(
      zone.getNucnetZone(),
      Libstatmech__Fermion__getRestMass( p_electron ),
      d_eta_F,
      zone.getElectronNeutrinoChemicalPotential()
    );

    zone.updateProperty(
      nnt::s_MUEKT,
      boost::lexical_cast<std::string>( d_mue_kT )
    );

    Libstatmech__Fermion__free( p_electron );

  }

  update_my_neutrino_rate_functions_data( zone );

}

//##############################################################################
// update_user_rate_functions_low_T_rates().
//##############################################################################

void
update_user_rate_functions_low_T_rates(
  Libnucnet__Reac * p_reac,
  const char * s_reac_xpath,
  const char * s_reac_file,
  double d_t9_cutoff,
  double d_t9_cutoff_factor
)
{

  //============================================================================
  // Get reaction view.
  //============================================================================

  Libnucnet__ReacView * p_view =
    Libnucnet__ReacView__new( p_reac, s_reac_xpath );

  //============================================================================
  // Get lab rates.
  //============================================================================

  Libnucnet__Reac * p_lab_reactions =
    Libnucnet__Reac__new_from_xml(
      s_reac_file, 
      ""
    );

  //============================================================================
  // Set transition to lab rates.
  //============================================================================

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list(
      Libnucnet__ReacView__getReac( p_view )
    );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    std::string s_key =
      Libnucnet__Reaction__getRateFunctionKey( reaction.getNucnetReaction() );

    if(
      s_key != SINGLE_RATE_STRING &&
      s_key != RATE_TABLE_STRING &&
      s_key != NON_SMOKER_STRING
    )
    {

      Libnucnet__Reaction * p_reaction =
	Libnucnet__Reac__getReactionByString(
	  p_lab_reactions,
	  Libnucnet__Reaction__getString( reaction.getNucnetReaction() )
	);

      if( p_reaction )
      {
	Libnucnet__Reaction__updateUserRateFunctionProperty(
	  reaction.getNucnetReaction(), 
	  nnt::s_LAB_RATE,
	  NULL,
	  NULL,
	  boost::lexical_cast<std::string>(
	    Libnucnet__Reaction__computeRate(
	      p_reaction,
	      d_t9_cutoff,
	      NULL
	    )
	  ).c_str()
	);
      }
      else
      {
	Libnucnet__Reaction__updateUserRateFunctionProperty(
	  reaction.getNucnetReaction(),
	  nnt::s_LAB_RATE,
	  NULL,
	  NULL,
	  "0."
	);
      }

      Libnucnet__Reaction__updateUserRateFunctionProperty(
	reaction.getNucnetReaction(),
	nnt::s_LAB_RATE_T9_CUTOFF,
	NULL,
	NULL,
	boost::lexical_cast<std::string>( d_t9_cutoff ).c_str()
      );

      Libnucnet__Reaction__updateUserRateFunctionProperty(
	reaction.getNucnetReaction(),
	nnt::s_LAB_RATE_T9_CUTOFF_FACTOR,
	NULL,
	NULL,
	boost::lexical_cast<std::string>( d_t9_cutoff_factor ).c_str()
      );

    }

  }

  //============================================================================
  // Clean up. Done.
  //============================================================================

  Libnucnet__ReacView__free( p_view );
  Libnucnet__Reac__free( p_lab_reactions );

}

} // namespace user
