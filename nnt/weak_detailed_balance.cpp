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
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//!
//! \file
//! \brief Routines to compute weak detailed balance.
//!
////////////////////////////////////////////////////////////////////////////////

#include "nnt/weak_detailed_balance.h"

namespace nnt
{

//##############################################################################
// compute_reverse_weak_rate_for_reaction().
//##############################################################################

double
compute_reverse_weak_rate_for_reaction(
  Libnucnet__Nuc *p_nuc,
  Libnucnet__Reaction *p_reaction,
  double d_forward,
  double d_t9,
  double d_rho,
  double d_muekT,
  double d_munuekT
)
{

  Libnucnet__Species * p_species;
  double d_reverse, d_exp = 0;

  if( d_munuekT == GSL_NEGINF ) return 0.; 

  if( GSL_SIGN( d_forward ) == GSL_SIGN( -d_forward ) ) return 0;

  reaction_element_list_t reactant_list =
    make_reaction_reactant_list( p_reaction );

  for(
    reaction_element_list_t::iterator iter = reactant_list.begin();
    iter != reactant_list.end();
    iter++
  )
  {

    if(
      Libnucnet__Reaction__Element__isNuclide(
	iter->getNucnetReactionElement()
      )
    )
    {

      p_species = 
	Libnucnet__Nuc__getSpeciesByName(
	  p_nuc,
	  Libnucnet__Reaction__Element__getName(
	    iter->getNucnetReactionElement()
	  )
	);

      d_exp +=
	(
	  log(
	    Libnucnet__Species__computeQuantumAbundance(
	      p_species,
	      d_t9,
	      d_rho
	    )
	  )
	  +
	  log( d_rho )
	  -
	  (
	    (int) Libnucnet__Species__getZ( p_species )
	  ) *
	  ( d_muekT - d_munuekT )
	  -
	  (
	    Libnucnet__Species__getMassExcess( p_species ) /
	    compute_kT_in_MeV( d_t9 )
	  )
	);

    }

  }

  reaction_element_list_t product_list =
    make_reaction_product_list( p_reaction );

  for(
    reaction_element_list_t::iterator iter = product_list.begin();
    iter != product_list.end();
    iter++
  )
  {

    if(
      Libnucnet__Reaction__Element__isNuclide(
	iter->getNucnetReactionElement()
      )
    )
    {

      p_species = 
	Libnucnet__Nuc__getSpeciesByName(
	  p_nuc,
	  Libnucnet__Reaction__Element__getName(
	    iter->getNucnetReactionElement()
	  )
	);

      d_exp -=
	(
	  log(
	    Libnucnet__Species__computeQuantumAbundance(
	      p_species,
	      d_t9,
	      d_rho
	    )
	  )
	  +
	  log( d_rho )
	  -
	  (
	    (int) Libnucnet__Species__getZ( p_species )
	  ) *
	  ( d_muekT - d_munuekT )
	  -
	  (
	    Libnucnet__Species__getMassExcess( p_species ) /
	    compute_kT_in_MeV( d_t9 )
	  )
	);

    } 

  }

  d_reverse =
    d_forward * exp( d_exp ) *
      Libnucnet__Reaction__getDuplicateProductFactor( p_reaction ) /
      Libnucnet__Reaction__getDuplicateReactantFactor( p_reaction );

  return d_reverse;

}

//##############################################################################
// Zone::setWeakDetailedBalance().
//##############################################################################

void
Zone::setWeakDetailedBalance( )
{

  Libnucnet__NetView * p_view, * p_evolution_view;
  Libstatmech__Fermion *p_electron;
  double d_mue_kT, d_forward, d_reverse;

  if(
    this->getProperty( s_MU_NUE_KT ) ==
    boost::lexical_cast<std::string>( GSL_NEGINF )
  )
    return;

  p_electron =
    Libstatmech__Fermion__new(
      s_ELECTRON,
      d_ELECTRON_MASS_IN_MEV,
      2,
      -1.
    );

  d_mue_kT =
    compute_electron_chemical_potential_kT(
      boost::lexical_cast<double>( this->getProperty( s_T9 ) ),
      boost::lexical_cast<double>( this->getProperty( s_RHO ) ),
      Libnucnet__Zone__computeZMoment( this->getNucnetZone(), 1 )
    );

  p_evolution_view = Libnucnet__Zone__getEvolutionNetView( getNucnetZone() );

  p_view = this->getNetView( "", s_WEAK_XPATH );

  reaction_list_t reaction_list =
    make_reaction_list(
      Libnucnet__Net__getReac( Libnucnet__NetView__getNet( p_view ) )
    );

  for(
    reaction_list_t::iterator iter = reaction_list.begin();
    iter != reaction_list.end();
    iter++
  )
  {

    if(
      Libnucnet__Reac__getReactionByString(
        Libnucnet__Net__getReac(
          Libnucnet__NetView__getNet( p_evolution_view )
        ),
        Libnucnet__Reaction__getString( iter->getNucnetReaction() )
      )
    )
    {       

      Libnucnet__Zone__getRatesForReaction(
	this->getNucnetZone(),
	iter->getNucnetReaction(),
	&d_forward,
	&d_reverse
      );

      d_reverse =
	compute_reverse_weak_rate_for_reaction(
	  Libnucnet__Net__getNuc(
	    Libnucnet__Zone__getNet( this->getNucnetZone() )
	  ),
	  iter->getNucnetReaction(),
	  d_forward,
	  boost::lexical_cast<double>( this->getProperty( s_T9 ) ),
	  boost::lexical_cast<double>( this->getProperty( s_RHO ) ),
	  d_mue_kT,
	  boost::lexical_cast<double>( this->getProperty( s_MU_NUE_KT ) )
	);

      if( !gsl_finite( d_reverse ) )
      {
	d_forward = 0.;
	d_reverse = 0.;
      }

      Libnucnet__Zone__updateRatesForReaction(
	this->getNucnetZone(),
	iter->getNucnetReaction(),
	d_forward,
	d_reverse
      );

    }

  }

  Libstatmech__Fermion__free( p_electron );

}

} // namespace nnt
