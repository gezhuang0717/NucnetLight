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
//! \brief Code to limit the evolution network.
//!
////////////////////////////////////////////////////////////////////////////////

#include "network_limiter.h"

/**
 * @brief A namespace for user-defined rate functions.
 */
namespace user
{

//##############################################################################
// check_reactants().
//##############################################################################

int
check_reactants(
  Libnucnet__Zone * p_zone,
  gsl_vector * p_abunds,
  Libnucnet__Reaction * p_reaction,
  double d_cutoff
)
{

  Libnucnet__Species * p_species;

  nnt::reaction_element_list_t reactant_list =
    nnt::make_reaction_nuclide_reactant_list( p_reaction );

  BOOST_FOREACH( nnt::ReactionElement reactant, reactant_list )
  {

    p_species =
      Libnucnet__Nuc__getSpeciesByName(
        Libnucnet__Net__getNuc( Libnucnet__Zone__getNet( p_zone ) ),
        Libnucnet__Reaction__Element__getName(
          reactant.getNucnetReactionElement()
        )
      );

    if( !p_species )
    {
      return 0;
    }

    if(
      gsl_vector_get(
        p_abunds,
        Libnucnet__Species__getIndex( p_species )
      ) < d_cutoff
    )
    {
      return 0;
    }

  }

  return 1;

}
  
//##############################################################################
// check_products().
//##############################################################################

int
check_products(
  Libnucnet__Zone * p_zone,
  gsl_vector * p_abunds,
  Libnucnet__Reaction * p_reaction,
  double d_cutoff
)
{

  Libnucnet__Species * p_species;

  nnt::reaction_element_list_t element_list =
    nnt::make_reaction_nuclide_product_list( p_reaction );

  BOOST_FOREACH( nnt::ReactionElement product, element_list )
  {

    p_species =
      Libnucnet__Nuc__getSpeciesByName(
        Libnucnet__Net__getNuc( Libnucnet__Zone__getNet( p_zone ) ),
        Libnucnet__Reaction__Element__getName(
          product.getNucnetReactionElement()
        )
      );

    if( !p_species )
    {
      return 0;
    }

    if(
      gsl_vector_get(
        p_abunds,
        Libnucnet__Species__getIndex( p_species )
      ) < d_cutoff
    )
    {
      return 0;
    }

  }

  return 1;

}
  
//##############################################################################
// limit_evolution_network().
//##############################################################################

/**
 * \brief Routine to limit reactions to those that connect species with
 *        abundance greater than a thresold value.
 *
 * \param zone The zone.
 * \param d_cutoff The cutoff threshold (optional--default = 1.e-25)
 *
*/

void
limit_evolution_network( nnt::Zone& zone, double d_cutoff )
{

  Libnucnet__NetView * p_base_view, * p_view;
  gsl_vector * p_abunds;

  if(
    zone.hasProperty( nnt::s_LOW_T_NETWORK_T9 ) &&
    zone.hasProperty( nnt::s_LOW_T_NETWORK_REAC_XPATH ) &&
    boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) <
    boost::lexical_cast<double>( zone.getProperty( nnt::s_LOW_T_NETWORK_T9 ) )
  )
  {
    p_base_view =
       zone.getNetView(
         "",
         zone.getProperty( nnt::s_LOW_T_NETWORK_REAC_XPATH).c_str()
       );
  }
  else
  {
    p_base_view = zone.getNetView( "", "" );
  }

  p_abunds = Libnucnet__Zone__getAbundances( zone.getNucnetZone() );

  boost::unordered_set<std::string> species_set;

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list(
      Libnucnet__Net__getReac( Libnucnet__NetView__getNet( p_base_view ) )
    );

  p_view = Libnucnet__NetView__copy( p_base_view );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    if(
      !check_reactants(
         zone.getNucnetZone(),
         p_abunds,
         reaction.getNucnetReaction(),
         d_cutoff
      )
      &&
      !check_products(
         zone.getNucnetZone(),
         p_abunds,
         reaction.getNucnetReaction(),
         d_cutoff
      )
    )
    {

      Libnucnet__NetView__removeReaction(
        p_view,
        reaction.getNucnetReaction()
      );

    }
    else
    {

      nnt::reaction_element_list_t reactant_list =
        nnt::make_reaction_nuclide_reactant_list(
          reaction.getNucnetReaction()
        );

      BOOST_FOREACH( nnt::ReactionElement element, reactant_list )
      {
        species_set.insert(
          Libnucnet__Reaction__Element__getName(
            element.getNucnetReactionElement()
          )  
        );
      }

      nnt::reaction_element_list_t product_list =
        nnt::make_reaction_nuclide_product_list(
          reaction.getNucnetReaction()
        );

      BOOST_FOREACH( nnt::ReactionElement element, product_list )
      {
        species_set.insert(
          Libnucnet__Reaction__Element__getName(
            element.getNucnetReactionElement()
          )  
        );
      }

    }

  }

  Libnucnet__Zone__updateNetView(
    zone.getNucnetZone(),
    EVOLUTION_NETWORK,
    NULL,
    NULL,
    p_view
  );

  nnt::species_list_t species_list =
    nnt::make_species_list(
      Libnucnet__Net__getNuc( Libnucnet__Zone__getNet( zone.getNucnetZone() ) )
    );

  BOOST_FOREACH( nnt::Species species, species_list )
  {

    if(
      !(
         species_set.find(
           Libnucnet__Species__getName( species.getNucnetSpecies() )
         ) != species_set.end()
      )
    )
    {

      if(
        gsl_vector_get(
          p_abunds,
          Libnucnet__Species__getIndex( species.getNucnetSpecies() )
        ) < d_cutoff
      )
        Libnucnet__Zone__updateSpeciesAbundance(
          zone.getNucnetZone(),
          species.getNucnetSpecies(),
          0.
        );
 
    }

  }
   
  gsl_vector_free( p_abunds );

}

void
limit_evolution_network( nnt::Zone& zone )
{

  limit_evolution_network( zone, 1.e-25 );

}

//##############################################################################
// zero_out_small_abundances().
//##############################################################################

/**
 * \brief Routine to zero out the abundances of species with
 *        abundance less than a thresold value.
 *
 * \param zone The zone.
 * \param d_threshold The cutoff threshold.
 *
*/

void
zero_out_small_abundances( nnt::Zone& zone, double d_threshold )
{

  nnt::species_list_t species_list =
    nnt::make_species_list(
      Libnucnet__Net__getNuc(
        Libnucnet__Zone__getNet( zone.getNucnetZone() )
      )
    );

  BOOST_FOREACH( nnt::Species species, species_list )
  {

    if(
      Libnucnet__Zone__getSpeciesAbundance(
        zone.getNucnetZone(),
        species.getNucnetSpecies()
      ) < d_threshold
    )
      Libnucnet__Zone__updateSpeciesAbundance(
        zone.getNucnetZone(),
        species.getNucnetSpecies(),
        0.
      );

  }

}

//##############################################################################
// zero_out_small_rates().
//##############################################################################

void
zero_out_small_rates( nnt::Zone& zone, double d_threshold )
{

/**
 * \brief Routine to zero out the rates for a reaction whose forward and
 *        reverse rates are less than a threshold value.
 *
 * \param zone The zone.
 * \param d_threshold The rate cutoff threshold.
 *
*/

  double d_forward, d_reverse;

  Libnucnet__NetView * p_view = zone.getNetView( EVOLUTION_NETWORK );
 
  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list(
      Libnucnet__Net__getReac(
        Libnucnet__NetView__getNet( p_view )
      )
    );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    Libnucnet__Zone__getRatesForReaction(
      zone.getNucnetZone(),
      reaction.getNucnetReaction(),
      &d_forward,
      &d_reverse
    );

    if( d_forward < d_threshold && d_reverse < d_threshold )
      Libnucnet__Zone__updateRatesForReaction(
        zone.getNucnetZone(),
        reaction.getNucnetReaction(),
        0.,
        0.
      );

  }

}
  
//##############################################################################
// get_isolated_species().
//##############################################################################

/**
 * \brief Routine to find species within a nuclear view that are isolated
 *        (no reaction in or out) according to the reaction view.
 *
 * \param p_net A Libnucnet__Net network.
 * \param s_nuc_xpath A string giving the nuclear view XPath expression.
 * \param s_reac_xpath A string giving the reaction view XPath expression.
 * \return A std::set containing strings giving the names of the isolated 
 *         species.
 *
*/

std::set<std::string>
get_isolated_species(
  Libnucnet__Net * p_net,
  std::string s_nuc_xpath,
  std::string s_reac_xpath
)
{

  Libnucnet__NetView * p_view;
  std::set<std::string> species_set;
  nnt::reaction_element_list_t element_list;

  p_view =
    Libnucnet__NetView__new(
      p_net,
      s_nuc_xpath.c_str(),
      s_reac_xpath.c_str()
    );

  Libnucnet__Nuc__setSpeciesCompareFunction(
    Libnucnet__Net__getNuc( p_net ),
    (Libnucnet__Species__compare_function) nnt::species_sort_function
  );

  Libnucnet__Nuc__sortSpecies(
    Libnucnet__Net__getNuc( p_net )
  );

  nnt::species_list_t species_list =
    nnt::make_species_list( Libnucnet__Net__getNuc( p_net ) );

  BOOST_FOREACH( nnt::Species species, species_list )
  {

    species_set.insert(
      Libnucnet__Species__getName( species.getNucnetSpecies() )
    );

  }

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list(
      Libnucnet__Net__getReac( Libnucnet__NetView__getNet( p_view ) )
    );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    element_list =
      nnt::make_reaction_nuclide_reactant_list( reaction.getNucnetReaction() );

    BOOST_FOREACH( nnt::ReactionElement element, element_list )
    {

      if(
        species_set.find(
          Libnucnet__Reaction__Element__getName(
            element.getNucnetReactionElement()
          )
        ) != species_set.end()
      )
      {
        species_set.erase(
          species_set.find(
            Libnucnet__Reaction__Element__getName(
              element.getNucnetReactionElement()
            )
          )
        );
      }

    }

    element_list =
      nnt::make_reaction_nuclide_product_list( reaction.getNucnetReaction() );

    BOOST_FOREACH( nnt::ReactionElement element, element_list )
    {

      if(
        species_set.find(
          Libnucnet__Reaction__Element__getName(
            element.getNucnetReactionElement()
          )
        ) != species_set.end()
      )
      {
        species_set.erase(
          species_set.find(
            Libnucnet__Reaction__Element__getName(
              element.getNucnetReactionElement()
            )
          )
        );
      }

    }

  }

  Libnucnet__NetView__free( p_view );

  return species_set;
      
}
  
} // namespace user
