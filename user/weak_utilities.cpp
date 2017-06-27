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
//! \file
//! \brief Code for useful utilities for weak reactions.
////////////////////////////////////////////////////////////////////////////////
  
#include "weak_utilities.h"

namespace user
{


//##############################################################################
// Global maps to store two-d weak quantities.
//##############################################################################

boost::unordered_map<std::string, nnt::TwoDWeakQuantity> my_rates;
boost::unordered_map<std::string, nnt::TwoDWeakQuantity> my_log10_fts;
boost::unordered_map<std::string, nnt::TwoDWeakQuantity> my_energy_loss;

//##############################################################################
// my_yedot().
//##############################################################################

double
my_yedot(
  double d_ye,
  void * p_data
)
{

  boost::tuple<double, double, double, double, double> t =
    compute_all_yedot(
      d_ye,
      p_data
    );

  return t.get<4>();

}
   
//##############################################################################
// compute_all_yedot().
//##############################################################################

boost::tuple<double, double, double, double, double>
compute_all_yedot(
  double d_ye,
  void * p_data
)
{

  double d_yedot_bm = 0;
  double d_yedot_bp = 0; 
  double d_yedot_ec = 0;
  double d_yedot_pc = 0;
  double d_yedot;

  nnt::Zone zone = *( nnt::Zone * ) p_data;

  //============================================================================
  // Update user rate functions.
  //============================================================================

  update_my_rate_functions_data( zone );

  //============================================================================
  // Update weak rate data.
  //============================================================================

  d_yedot_bm =
    compute_total_flow(
      zone,
      nnt::s_FORWARD_FLOW,
      zone.getNetView( "", nnt::s_BETA_MINUS_XPATH )
    ) -
    compute_total_flow(
      zone,
      nnt::s_REVERSE_FLOW,
      zone.getNetView( "", nnt::s_BETA_MINUS_XPATH )
    );

  d_yedot_bp =
    compute_total_flow(
      zone,
      nnt::s_FORWARD_FLOW,
      zone.getNetView( "", nnt::s_BETA_PLUS_XPATH )
    ) -
    compute_total_flow(
      zone,
      nnt::s_REVERSE_FLOW,
      zone.getNetView( "", nnt::s_BETA_PLUS_XPATH )
    );

  d_yedot_ec =
    compute_total_flow(
      zone,
      nnt::s_FORWARD_FLOW,
      zone.getNetView( "", nnt::s_ELECTRON_CAPTURE_XPATH )
    ) -
    compute_total_flow(
      zone,
      nnt::s_REVERSE_FLOW,
      zone.getNetView( "", nnt::s_ELECTRON_CAPTURE_XPATH )
    );

  d_yedot_pc =
    compute_total_flow(
      zone,
      nnt::s_FORWARD_FLOW,
      zone.getNetView( "", nnt::s_POSITRON_CAPTURE_XPATH )
    ) -
    compute_total_flow(
      zone,
      nnt::s_REVERSE_FLOW,
      zone.getNetView( "", nnt::s_POSITRON_CAPTURE_XPATH )
    );

  d_yedot = d_yedot_bm - d_yedot_bp - d_yedot_ec + d_yedot_pc;

  return
    boost::make_tuple(
      d_yedot_bm,
      d_yedot_bp,
      d_yedot_ec,
      d_yedot_pc,
      d_yedot
    );
    
}

//##############################################################################
// compute_dyedot_dye().
//##############################################################################

double
compute_dyedot_dye(
  nnt::Zone& zone,
  double d_variation
)
{

  gsl_function F;
  double d_result, d_abserr;
     
  F.function = &my_yedot_root;
  F.params = &zone;


  gsl_deriv_central(
    &F,
    Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 ),
    d_variation,
    &d_result,
    &d_abserr
  );

  return d_result;

}

//##############################################################################
// set_weak_views_in_zones().
//##############################################################################

void
set_weak_views_in_zones(
  Libnucnet * p_nucnet
)
{

  Libnucnet__NetView * p_view_bm, * p_view_bp, * p_view_ec, * p_view_pc;

  p_view_bm =
    Libnucnet__NetView__new(
      Libnucnet__getNet( p_nucnet ),
      "",
      nnt::s_BETA_MINUS_XPATH
    );

  p_view_bp =
    Libnucnet__NetView__new(
      Libnucnet__getNet( p_nucnet ),
      "",
      nnt::s_BETA_PLUS_XPATH
    );

  p_view_ec =
    Libnucnet__NetView__new(
      Libnucnet__getNet( p_nucnet ),
      "",
      nnt::s_ELECTRON_CAPTURE_XPATH
    );

  p_view_pc =
    Libnucnet__NetView__new(
      Libnucnet__getNet( p_nucnet ),
      "",
      nnt::s_POSITRON_CAPTURE_XPATH
    );

  nnt::zone_list_t zone_list = nnt::make_zone_list( p_nucnet );

  for(
    nnt::zone_list_t::iterator it = zone_list.begin();
    it != zone_list.end();
    it++
  )
  {
  
    Libnucnet__Zone__updateNetView(
      it->getNucnetZone(),
      "",
      nnt::s_BETA_MINUS_XPATH,
      NULL,
      Libnucnet__NetView__copy( p_view_bm )
    );

    Libnucnet__Zone__updateNetView(
      it->getNucnetZone(),
      "",
      nnt::s_BETA_PLUS_XPATH,
      NULL,
      Libnucnet__NetView__copy( p_view_bp )
    );

    Libnucnet__Zone__updateNetView(
      it->getNucnetZone(),
      "",
      nnt::s_ELECTRON_CAPTURE_XPATH,
      NULL,
      Libnucnet__NetView__copy( p_view_ec )
    );

    Libnucnet__Zone__updateNetView(
      it->getNucnetZone(),
      "",
      nnt::s_POSITRON_CAPTURE_XPATH,
      NULL,
      Libnucnet__NetView__copy( p_view_pc )
    );

  }

  Libnucnet__NetView__free( p_view_bm );
  Libnucnet__NetView__free( p_view_bp );
  Libnucnet__NetView__free( p_view_ec );
  Libnucnet__NetView__free( p_view_pc );

}

//##############################################################################
// my_yedot_root().
//##############################################################################

double
my_yedot_root(
  double d_ye,
  void *p_data
)
{

  nnt::Zone * p_nnt_zone;

  //============================================================================
  // Ye minimum.
  //============================================================================

  if( d_ye < 1.e-300 ) d_ye = 1.e-300;

  //============================================================================
  // Ye maximum.
  //============================================================================

  if( d_ye > 1. ) d_ye = 1.;

  //============================================================================
  // Data.
  //============================================================================

  p_nnt_zone = ( nnt::Zone * ) p_data;

  //============================================================================
  // Update abundances with Equilibrium.
  //============================================================================

  p_nnt_zone->updateProperty(
    nnt::s_YE,
    boost::lexical_cast<std::string>( d_ye )
  );

  nnt::set_zone_abundances_to_equilibrium( *p_nnt_zone );

  //============================================================================
  // Return Yedot.
  //============================================================================

  return
    my_yedot( d_ye, p_nnt_zone );

}

//##############################################################################
// compute_lambda().
//##############################################################################

double
compute_lambda( nnt::Zone zone )
{

  return
    (
      boost::lexical_cast<double>( zone.getProperty( nnt::s_MUPKT ) ) +
      (
	Libnucnet__Species__getMassExcess(
	  Libnucnet__Nuc__getSpeciesByName(
	    Libnucnet__Net__getNuc(
	      Libnucnet__Zone__getNet( zone.getNucnetZone() )
	    ),
	    "h1"
	  )
	)
	/
	nnt::compute_kT_in_MeV(
	  boost::lexical_cast<double>(
            zone.getProperty( nnt::s_T9 )
	  )
	)
      )
      +
      boost::lexical_cast<double>(
        zone.getProperty( nnt::s_MUEKT )
      ) -
      boost::lexical_cast<double>(
        zone.getProperty( nnt::s_MUNKT )
      ) -
      (
	Libnucnet__Species__getMassExcess(
	  Libnucnet__Nuc__getSpeciesByName(
	    Libnucnet__Net__getNuc(
	      Libnucnet__Zone__getNet( zone.getNucnetZone() )
	    ),
	    "n"
	  )
	)
	/
	nnt::compute_kT_in_MeV(
	  boost::lexical_cast<double>(
            zone.getProperty( nnt::s_T9 )
	  )
	)
      )
    )
    /
    ( 
      boost::lexical_cast<double>(
        zone.getProperty( S_DYEDOT_DYE_NUCLEON )
      )
      -
      boost::lexical_cast<double>(
        zone.getProperty( S_DYEDOT_DYE )
      )
    );

}

//##############################################################################
// compute_approximate_neutrino_entropy_loss_rate().
//##############################################################################

double
compute_approximate_neutrino_entropy_loss_rate(
  nnt::Zone zone
)
{

  Libnucnet__NetView * p_view;
  Libstatmech__Fermion * p_electron;
  double d_eta_F, d_result = 0, d_entropy_flow;

  p_view =
    zone.getNetView(
      "",
      "[ \
         (product = 'electron' and product = 'anti-neutrino_e') or \
         (reactant = 'electron' and product = 'neutrino_e') or \
         (product = 'positron' and product = 'neutrino_e') or \
         (reactant = 'positron' and product = 'anti-neutrino_e') \
       ]"
    );

  p_electron =
    Libstatmech__Fermion__new(
      nnt::s_ELECTRON,
      nnt::d_ELECTRON_MASS_IN_MEV,
      2,
      -1.
    );

  d_eta_F =
    Libstatmech__Fermion__computeChemicalPotential(
      p_electron,
      boost::lexical_cast<double>(
        zone.getProperty( nnt::s_T9 )
      ) * GSL_CONST_NUM_GIGA,
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) *
        Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 ) *
        GSL_CONST_NUM_AVOGADRO,
      NULL,
      NULL
    )
    +
    Libstatmech__Fermion__getRestMass( p_electron ) /
      nnt::compute_kT_in_MeV(
        boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) )
      );

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list(
      Libnucnet__Net__getReac(
        Libnucnet__NetView__getNet( p_view )
      )
    );

  for(
    nnt::reaction_list_t::iterator it = reaction_list.begin();
    it != reaction_list.end();
    it++
  )
  {

    nnt::reaction_element_list_t reactant_list =
      nnt::make_reaction_nuclide_reactant_list( it->getNucnetReaction() );

    nnt::reaction_element_list_t product_list =
      nnt::make_reaction_nuclide_product_list( it->getNucnetReaction() );

    if(
      reactant_list.size() == 1 &&
      product_list.size() == 1
    )
    {

      d_entropy_flow =
        aa522a25__compute_reaction_neutrino_energy_loss_rate(
	  it->getNucnetReaction(),
	  Libnucnet__Zone__getNet( zone.getNucnetZone() ),
	  Libstatmech__Fermion__getRestMass( p_electron ),
	  boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ),
	  d_eta_F,
          zone.getElectronNeutrinoChemicalPotential()
        );

      for(
	nnt::reaction_element_list_t::iterator
	  it_reac = reactant_list.begin();
	it_reac != reactant_list.end();
	it_reac++
      )
      {

	d_entropy_flow *=
	  Libnucnet__Zone__getSpeciesAbundance(
	    zone.getNucnetZone(),
	    Libnucnet__Nuc__getSpeciesByName(
	      Libnucnet__Net__getNuc(
		Libnucnet__Zone__getNet( zone.getNucnetZone() )
	      ),
	      Libnucnet__Reaction__Element__getName(
		it_reac->getNucnetReactionElement()
	      )
	    )
	  );

      }

      d_result += d_entropy_flow;

    }

  }

  Libstatmech__Fermion__free( p_electron );

  return d_result;

}  

//##############################################################################
// compute_neutrino_entropy_loss_rate().
//##############################################################################

double 
compute_neutrino_entropy_loss_rate(
  nnt::Zone zone
)
{

  //============================================================================
  // Only compute loss if neutrinos escape.
  //============================================================================

  if( zone.getProperty( nnt::s_MU_NUE_KT ) != "-inf" )
    return 0;
  else
    return
      compute_neutrino_energy_loss_rate( zone ) /
      nnt::compute_kT_in_MeV(
        boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) )
      );

}

//##############################################################################
// compute_neutrino_energy_loss_rate().
//##############################################################################

double
compute_neutrino_energy_loss_rate(
  nnt::Zone zone
)
{

  Libnucnet__NetView * p_view;
  Libstatmech__Fermion * p_electron;
  double d_eta_F, d_energy_flow, d_result = 0.;

  //============================================================================
  // Get reaction view that emit neutrino_e and anti-neutrino_e.
  //============================================================================

  p_view =
    zone.getNetView(
      "",
      "[product = 'neutrino_e' or product = 'anti-neutrino_e']" 
    ); 

  //============================================================================
  // Set reaction compare function. 
  //============================================================================

  Libnucnet__Reac__setReactionCompareFunction(
    Libnucnet__Net__getReac( Libnucnet__NetView__getNet( p_view ) ),
    (Libnucnet__Reaction__compare_function) nnt::compare_reactions_by_string
  );

  //============================================================================
  // Set parameters. 
  //============================================================================

  p_electron =
    Libstatmech__Fermion__new(
      nnt::s_ELECTRON,
      nnt::d_ELECTRON_MASS_IN_MEV,
      2,
      -1.
    );

  zone.updateProperty(
    nnt::s_YE,
    boost::lexical_cast<std::string>(
      Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 )
    )
  );

  d_eta_F =
    Libstatmech__Fermion__computeChemicalPotential(
      p_electron,
      boost::lexical_cast<double>(
        zone.getProperty( nnt::s_T9 ) 
      ) * GSL_CONST_NUM_GIGA,
      boost::lexical_cast<double>(
        zone.getProperty( nnt::s_RHO ) 
      ) *
        boost::lexical_cast<double>( zone.getProperty( nnt::s_YE ) ) *
        GSL_CONST_NUM_AVOGADRO,
      NULL,
      NULL
    )
    +
    Libstatmech__Fermion__getRestMass( p_electron ) /
      nnt::compute_kT_in_MeV(
        boost::lexical_cast<double>(
          zone.getProperty( nnt::s_T9 )
        )
      );

  //============================================================================
  // Iterate reactions.
  //============================================================================
  
  nnt::reaction_list_t reaction_list = 
    nnt::make_reaction_list( 
      Libnucnet__Net__getReac( 
        Libnucnet__NetView__getNet( p_view )
      )
    );

  for(
    nnt::reaction_list_t::iterator it = reaction_list.begin();
    it != reaction_list.end();
    it++
  ) 
  {

    nnt::reaction_element_list_t reactant_list = 
      nnt::make_reaction_nuclide_reactant_list( it->getNucnetReaction() );

    nnt::reaction_element_list_t product_list = 
      nnt::make_reaction_nuclide_product_list( it->getNucnetReaction() );

    if( 
      reactant_list.size() == 1 &&
      product_list.size() == 1
    )
    {

      d_energy_flow = 
        compute_reaction_neutrino_energy_loss_rate(
          it->getNucnetReaction(),
          Libnucnet__Zone__getNet( zone.getNucnetZone() ),
          Libstatmech__Fermion__getRestMass( p_electron ),
          boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ),
          boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) *
            boost::lexical_cast<double>( zone.getProperty( nnt::s_YE ) ),
          zone.getElectronNeutrinoChemicalPotential(),
          d_eta_F
        );

      for(
        nnt::reaction_element_list_t::iterator it_reac = reactant_list.begin();
        it_reac != reactant_list.end();
        it_reac++
      )
      {

        d_energy_flow *=
          Libnucnet__Zone__getSpeciesAbundance(
            zone.getNucnetZone(),
            Libnucnet__Nuc__getSpeciesByName(
              Libnucnet__Net__getNuc(
                Libnucnet__Zone__getNet( zone.getNucnetZone() )
              ),
              Libnucnet__Reaction__Element__getName(
                it_reac->getNucnetReactionElement()
              )
            )
          );

      } 

      d_result += d_energy_flow;

    }

  }  

  Libstatmech__Fermion__free( p_electron );

  return d_result;

}

//##############################################################################
// compute_reaction_neutrino_energy_loss_rate().
//##############################################################################

double
compute_reaction_neutrino_energy_loss_rate(
  Libnucnet__Reaction *p_reaction,
  Libnucnet__Net *p_my_net,
  double d_electron_mass,
  double d_t9,
  double d_rhoe,
  double d_mu_nue_kT,
  double d_eta_F
)
{

  double d_average_energy;
  double d_rate = 0;
  double d_energy_loss_rate;
  char s_property[32];

  if( !Libnucnet__Net__isValidReaction( p_my_net, p_reaction ) )
    return 0.;

  Libnucnet__Reaction__iterateProducts(
    p_reaction,
    (Libnucnet__Reaction__Element__iterateFunction) set_neutrino_type,
    s_property
  );

  if(
    strcmp(
      Libnucnet__Reaction__getSource( p_reaction ),
      "fullerdat.txt"
    ) == 0
  )
  {

    boost::unordered_map<std::string, nnt::TwoDWeakQuantity>::iterator
      it_average_e =
        my_energy_loss.find(
          Libnucnet__Reaction__getString( p_reaction )
        );

    if( it_average_e != my_energy_loss.end() )
    {
      d_average_energy =
        it_average_e->second.computeValue( d_t9, d_rhoe ).first;
    }
    else
    {
      nnt::TwoDWeakQuantity my_average_energy( p_reaction, s_property );
      d_average_energy =
        my_average_energy.computeValue( d_t9, d_rhoe ).first;
    }

    if(
      strcmp(
        Libnucnet__Reaction__getRateFunctionKey( p_reaction ), 
        "two-d weak rates log10 ft"
      ) == 0
    )
    {
      d_rate = 
        compute_weak_rate_from_log10_ft(
          p_reaction,
          p_my_net,
          d_electron_mass,
          d_t9,
          d_rhoe,
          d_eta_F,
          d_mu_nue_kT
        );
    }
    else if(
      strcmp(
        Libnucnet__Reaction__getRateFunctionKey( p_reaction ), 
        "two-d weak rates"
      ) == 0
    )
    {
      d_rate = 
        compute_two_d_weak_rate(
          p_reaction,
          d_t9, 
          &d_rhoe
        );
    }

    return d_average_energy * d_rate;

  }
  else if(
    strcmp(
      Libnucnet__Reaction__getSource( p_reaction ),
      "aa522a25"
    ) == 0
  )
  {

    d_energy_loss_rate =
      aa522a25__compute_reaction_neutrino_energy_loss_rate(
        p_reaction,
        p_my_net,
        d_electron_mass, 
        d_t9,
        d_eta_F,
        d_mu_nue_kT
      );

    return d_energy_loss_rate;

  }
  else
  {
 
    fprintf( stderr, "No neutrino energy info in the reaction.\n" );
    return 0.;

  } 
 
}
  
//##############################################################################
// compute_approximate_neutrino_entropy_loss_rate().
//##############################################################################

double
compute_approximate_neutrino_entropy_loss_rate(
  nnt::Zone &zone
)
{

  Libnucnet__NetView * p_view;
  double d_result = 0.;
  std::pair< double, double > flows;
  double d_kT, d_eta_F;
  Libstatmech__Fermion * p_electron;

  //============================================================================
  // Check that neutrinos lost.
  //============================================================================

  if( zone.getProperty( nnt::s_MU_NUE_KT ) != "-inf" ) return 0;

  //============================================================================
  // Get reaction view that emit neutrino_e and anti-neutrino_e.
  //============================================================================

  p_view =
    zone.getNetView(
      "",
      "[product = 'neutrino_e' or product = 'anti-neutrino_e']" 
    ); 

  //============================================================================
  // Get kT.
  //============================================================================

  d_kT =
    nnt::compute_kT_in_MeV(
      boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) )
    );
  
  //============================================================================
  // Get electron chemical potential eta_F.
  //============================================================================

  p_electron =
    Libstatmech__Fermion__new(
      nnt::s_ELECTRON,
      nnt::d_ELECTRON_MASS_IN_MEV,
      2,
      -1.
    );

  d_eta_F =
    Libstatmech__Fermion__computeChemicalPotential(
      p_electron,
      boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) *
        GSL_CONST_NUM_GIGA,
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) *
        Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 ) *
          GSL_CONST_NUM_AVOGADRO,
      NULL,
      NULL
    )
    +
    Libstatmech__Fermion__getRestMass( p_electron ) / d_kT;

  //============================================================================
  // Iterate reactions.
  //============================================================================
  
  nnt::reaction_list_t reaction_list = 
    nnt::make_reaction_list( 
      Libnucnet__Net__getReac( 
        Libnucnet__NetView__getNet( p_view )
      )
    );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    flows =
      compute_flows_for_reaction(
        zone,
        reaction.getNucnetReaction()
      );

    nnt::reaction_element_list_t reactant_list = 
      nnt::make_reaction_nuclide_reactant_list( reaction.getNucnetReaction() );

    nnt::reaction_element_list_t product_list = 
      nnt::make_reaction_nuclide_product_list( reaction.getNucnetReaction() );

    if( 
      reactant_list.size() == 1 &&
      product_list.size() == 1
    )
    {

      d_result +=
        flows.first *
        aa522a25__compute_reaction_average_neutrino_energy(
          reaction.getNucnetReaction(),
          Libnucnet__Zone__getNet( zone.getNucnetZone() ),
          Libstatmech__Fermion__getRestMass( p_electron ),
          boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ),
          d_eta_F,
          zone.getElectronNeutrinoChemicalPotential()
        ) / d_kT;

    }

  }  

  Libstatmech__Fermion__free( p_electron );

  return d_result;

}

//##############################################################################
// compute_simple_neutrino_entropy_loss_rate().
//##############################################################################

double
compute_simple_neutrino_entropy_loss_rate(
  nnt::Zone zone,
  double d_neutrino_factor
)
{

  Libnucnet__NetView * p_view;
  double d_result = 0.;
  std::pair< double, double > flows;

  //============================================================================
  // Check that neutrinos lost.
  //============================================================================

  if( zone.getProperty( nnt::s_MU_NUE_KT ) != "-inf" ) return 0;

  //============================================================================
  // Get reaction view that emit neutrino_e and anti-neutrino_e.
  //============================================================================

  p_view =
    zone.getNetView(
      "",
      "[product = 'neutrino_e' or product = 'anti-neutrino_e']" 
    ); 

  //============================================================================
  // Iterate reactions.
  //============================================================================
  
  nnt::reaction_list_t reaction_list = 
    nnt::make_reaction_list( 
      Libnucnet__Net__getReac( 
        Libnucnet__NetView__getNet( p_view )
      )
    );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    flows =
      compute_flows_for_reaction(
        zone,
        reaction.getNucnetReaction()
      );

    nnt::reaction_element_list_t reactant_list = 
      nnt::make_reaction_nuclide_reactant_list( reaction.getNucnetReaction() );

    nnt::reaction_element_list_t product_list = 
      nnt::make_reaction_nuclide_product_list( reaction.getNucnetReaction() );

    if( 
      reactant_list.size() == 1 &&
      product_list.size() == 1
    )
    {

      d_result +=
        d_neutrino_factor *
        flows.first *
        nnt::compute_reaction_nuclear_Qvalue(
          Libnucnet__NetView__getNet( p_view ),
          reaction.getNucnetReaction(),
          nnt::d_ELECTRON_MASS_IN_MEV
        ) /
        nnt::compute_kT_in_MeV(
          boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) )
        );

    }

  }  

  return d_result;

}

//##############################################################################
// set_neutrino_type().
//##############################################################################

int
set_neutrino_type(
  Libnucnet__Reaction__Element *p_element,
  char *s_property
)
{

  if( 
    strcmp(
      Libnucnet__Reaction__Element__getName( p_element ),
      "neutrino_e"
    ) == 0
  ) 
  {
    strcpy( s_property, nnt::s_AVERAGE_ENERGY_NU_E );
  }
  else if( 
    strcmp(
      Libnucnet__Reaction__Element__getName( p_element ),
      "anti-neutrino_e"
    ) == 0
  ) 
  {
    strcpy( s_property, nnt::s_AVERAGE_ENERGY_NUBAR_E );
  }

  return 1;

}

//############################################################################
// compute_weak_rate_from_log10_ft().
//############################################################################

double
compute_weak_rate_from_log10_ft(
  Libnucnet__Reaction *p_reaction,
  Libnucnet__Net * p_net,
  double d_electron_mass,
  double d_t9,
  double d_rhoe,
  double d_etaF,
  double d_mu_nue_kT
)
{

  double d_result;

  std::pair<double,double> my_pair;

  boost::unordered_map<std::string,nnt::TwoDWeakQuantity>::iterator it =
    my_log10_fts.find(
      Libnucnet__Reaction__getString( p_reaction )
    );

  if( it != my_log10_fts.end() )
  {
    my_pair = it->second.computeValue( d_t9, d_rhoe );
  }
  else
  {
    nnt::TwoDWeakQuantity my_log10_ft( p_reaction, nnt::s_LOG10_FT );
    my_pair = my_log10_ft.computeValue( d_t9, d_rhoe );
  }

  if( my_pair.second < 2. )
  {

    d_result =
      M_LN2 *
      nnt::ffnIV__compute_Ie( 
        p_reaction, 
        p_net,
        d_electron_mass,
        d_t9, 
        d_etaF,
        d_mu_nue_kT
      ) /
      pow(
        10.,
        my_pair.first
      ); 

      correct_for_weak_lab_rate( p_reaction, d_t9, d_result );

  }
  else
  {

    d_result =
      compute_two_d_weak_rate(
        p_reaction,
        d_t9,
        &d_rhoe
      );

  }

  return d_result;
 
}

//##############################################################################
// set_two_d_weak_rates_hashes().
//##############################################################################

void
set_two_d_weak_rates_hashes( Libnucnet__Reac * p_reac )
{

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
        nnt::s_TWO_D_WEAK_RATES
      ) == 0
    )
    {
      my_rates.insert(
        std::make_pair(
          Libnucnet__Reaction__getString( it->getNucnetReaction() ),
          nnt::TwoDWeakQuantity( it->getNucnetReaction(), nnt::s_LOG10_RATE )
        )
      );
    }
    else if(
      strcmp(
        Libnucnet__Reaction__getRateFunctionKey( it->getNucnetReaction() ),
        nnt::s_TWO_D_WEAK_RATES_LOG10_FT
      ) == 0
    )
    {
      my_log10_fts.insert(
        std::make_pair(
          Libnucnet__Reaction__getString( it->getNucnetReaction() ),
          nnt::TwoDWeakQuantity( it->getNucnetReaction(), nnt::s_LOG10_FT )
        )
      );
      my_rates.insert(
        std::make_pair(
          Libnucnet__Reaction__getString( it->getNucnetReaction() ),
          nnt::TwoDWeakQuantity( it->getNucnetReaction(), nnt::s_LOG10_RATE )
        )
      );
    }

  }

}
    
//##############################################################################
// set_two_d_weak_energy_loss_hash().
//##############################################################################

void
set_two_d_weak_energy_loss_hash( Libnucnet__Reac * p_reac )
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
        nnt::s_TWO_D_WEAK_RATES
      ) == 0 ||
      strcmp(
        Libnucnet__Reaction__getRateFunctionKey( it->getNucnetReaction() ),
        nnt::s_TWO_D_WEAK_RATES_LOG10_FT
      ) == 0
    )
    {

      i_size = 0;

      Libnucnet__Reaction__iterateUserRateFunctionProperties(
        it->getNucnetReaction(),
        nnt::s_AVERAGE_ENERGY_NU_E,
        NULL,
        NULL,
        (Libnucnet__Reaction__user_rate_property_iterate_function)
           nnt::get_property_array_size,
        &i_size
      );

      if( i_size > 0 )
      {
        my_energy_loss.insert(
          std::make_pair(
            Libnucnet__Reaction__getString( it->getNucnetReaction() ),
            nnt::TwoDWeakQuantity(
               it->getNucnetReaction(),
               nnt::s_AVERAGE_ENERGY_NU_E
            )
          )
        );
      }

      i_size = 0;

      Libnucnet__Reaction__iterateUserRateFunctionProperties(
        it->getNucnetReaction(),
        nnt::s_AVERAGE_ENERGY_NUBAR_E,
        NULL,
        NULL,
        (Libnucnet__Reaction__user_rate_property_iterate_function)
           nnt::get_property_array_size,
        &i_size
      );

      if( i_size > 0 )
      {
        my_energy_loss.insert(
          std::make_pair(
            Libnucnet__Reaction__getString( it->getNucnetReaction() ),
            nnt::TwoDWeakQuantity(
               it->getNucnetReaction(),
               nnt::s_AVERAGE_ENERGY_NUBAR_E
            )
          )
        );
      }

    }

  }

}
    
//##############################################################################
// my_log10_ft_rate_function().
//##############################################################################

double
my_log10_ft_rate_function(
  Libnucnet__Reaction * p_reaction,
  double d_t9,
  void * p_data
)
{

  typedef struct {
    Libnucnet__Net *pNet;
    double dElectronMass;
    double dRhoe;
    double dEtaF;
    double dMuNuekT;
  } work;

  work * p_work = ( work * ) p_data;

  if( !p_data )
  {
    std::cerr << "Invalid data structure in my_log10_ft_rate_function."
              << std::endl;
    exit( EXIT_FAILURE );
  }

  return
    compute_weak_rate_from_log10_ft(
      p_reaction,
      p_work->pNet,
      p_work->dElectronMass,
      d_t9,
      p_work->dRhoe,
      p_work->dEtaF,
      p_work->dMuNuekT
    );

}

//##############################################################################
// my_approximate_weak_rate_function().
//##############################################################################

double
my_approximate_weak_rate_function(
  Libnucnet__Reaction * p_reaction,
  double d_t9,
  void * p_data
)
{
  
  typedef struct {
    Libnucnet__Net *pNet;
    double dElectronMassMeV;
    double dEtaF;
    double dMuNuekT;
  } work;

  double d_result;

  work * p_work = ( work * ) p_data;

  if( !p_data )
  {
    std::cerr << "Invalid data structure in my_approximate_weak_rate_function."
              << std::endl;
    exit( EXIT_FAILURE );
  }


  d_result =
    aa522a25__compute_rate(
      p_reaction,
      p_work->pNet,
      p_work->dElectronMassMeV,
      d_t9,
      p_work->dEtaF,
      p_work->dMuNuekT
    );

  correct_for_weak_lab_rate( p_reaction, d_t9, d_result );

  return d_result;

}

//##############################################################################
// compute_two_d_weak_rate().
//##############################################################################

double
compute_two_d_weak_rate(
  Libnucnet__Reaction * p_reaction,
  double d_t9,
  double *p_rhoe
)
{

  double d_result;

  if( !p_rhoe )
  {
    std::cerr << "Rhoe not set in two-d weak rate function." << std::endl;
    exit( EXIT_FAILURE );
  }

  boost::unordered_map<std::string,nnt::TwoDWeakQuantity>::iterator
    it = my_rates.find( Libnucnet__Reaction__getString( p_reaction ) );

  if( it != my_rates.end() )
  {
    d_result = it->second.computeValue( d_t9, *p_rhoe ).first;
  }
  else
  {
    nnt::TwoDWeakQuantity my_log10_rate( p_reaction, nnt::s_LOG10_RATE );
    d_result = my_log10_rate.computeValue( d_t9, *p_rhoe ).first;
  }

  if( d_result < -50. )
    d_result = 0.;
  else
    d_result = pow( 10., d_result );

  correct_for_weak_lab_rate( p_reaction, d_t9, d_result );

  return d_result;

}

//##############################################################################
// correct_for_weak_lab_rate().
//##############################################################################

void
correct_for_weak_lab_rate(
  Libnucnet__Reaction * p_reaction,
  double d_t9,
  double &d_rate
)
{

  const char * s_lab_rate, * s_t9_cutoff, * s_t9_cutoff_factor;
  double d_t9_cutoff, d_t9_cutoff_factor, d_f;

  s_lab_rate =
    Libnucnet__Reaction__getUserRateFunctionProperty(
      p_reaction,
      nnt::s_LAB_RATE,
      NULL,
      NULL
    );

  if( s_lab_rate )
  {

    s_t9_cutoff =
      Libnucnet__Reaction__getUserRateFunctionProperty(
        p_reaction,
        nnt::s_LAB_RATE_T9_CUTOFF,
        NULL,
        NULL
      );

    if( !s_t9_cutoff )
      std::cerr << "No cutoff temperature for lab rate." << std::endl;

    d_t9_cutoff = boost::lexical_cast<double>( s_t9_cutoff );

    s_t9_cutoff_factor =
      Libnucnet__Reaction__getUserRateFunctionProperty(
        p_reaction,
        nnt::s_LAB_RATE_T9_CUTOFF_FACTOR,
        NULL,
        NULL
      );

    if( !s_t9_cutoff_factor )
      std::cerr << "No cutoff factor for lab rate." << std::endl;

    d_t9_cutoff_factor = boost::lexical_cast<double>( s_t9_cutoff_factor );

    if( d_t9_cutoff > 0 )
    {
      d_f =
        boost::math::erfc(
          (
            log10( d_t9 / d_t9_cutoff ) / d_t9_cutoff_factor
          )
        ) / 2.;
    }
    else
    {
      d_f = 0.;
    }

    d_rate =
      ( 1. - d_f ) * d_rate +
	  d_f * boost::lexical_cast<double>( s_lab_rate );

  }

}

} // namespace user
