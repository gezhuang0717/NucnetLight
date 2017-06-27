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
//! \brief Code for routines related to reaction flows.
//!
////////////////////////////////////////////////////////////////////////////////

#include "user/flow_utilities.h"

namespace user
{

//############################################################################
// compute_flows_for_reaction().
//##########################################################################//

/**
 * \brief Compute the forward and reverse flows for a reaction in a zone.
 *
 * \param zone A Nucnet Tools zone.
 * \param p_reaction A pointer to the Libnucnet__Reaction.
 * \return A pair in which the first element is the forward flow and the
 *         second is the reverse.
 */

std::pair<double,double>
compute_flows_for_reaction(
  nnt::Zone &zone,
  Libnucnet__Reaction *p_reaction
)
{

  double d_forward_rate, d_reverse_rate;
  double d_f, d_r;
  size_t i_elements;

  Libnucnet__Net__computeRatesForReaction(
    Libnucnet__Zone__getNet( zone.getNucnetZone() ),
    p_reaction,
    boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ),
    boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ),
    Libnucnet__Zone__getDataForUserRateFunction(
      zone.getNucnetZone(),
      Libnucnet__Reaction__getRateFunctionKey( p_reaction )
    ),
    &d_forward_rate,
    &d_reverse_rate
  );

  if( Libnucnet__Reaction__isWeak( p_reaction ) )
  {

    d_reverse_rate =
      nnt::compute_reverse_weak_rate_for_reaction(
	Libnucnet__Net__getNuc(
          Libnucnet__Zone__getNet( zone.getNucnetZone() )
        ),
	p_reaction,
	d_forward_rate,
	boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ),
	boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ),
	boost::lexical_cast<double>( zone.getProperty( nnt::s_MUEKT ) ),
        zone.getElectronNeutrinoChemicalPotential()
      );

  }

  //==========================================================================
  // Modify rates for reaction.
  //==========================================================================

  modify_rates_for_reaction(
    zone,
    p_reaction,
    d_forward_rate,
    d_reverse_rate
  );

  //==========================================================================
  // Iterate the reactants to get the forward flow.
  //==========================================================================

  nnt::reaction_element_list_t reactant_list =
    nnt::make_reaction_nuclide_reactant_list( p_reaction );

  i_elements = reactant_list.size();

  d_f =
    d_forward_rate *
    pow(
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ),
      (double) i_elements - 1.
    ) /
      Libnucnet__Reaction__getDuplicateReactantFactor( p_reaction );

  BOOST_FOREACH( nnt::ReactionElement reactant, reactant_list )
  {

    d_f *=
      Libnucnet__Zone__getSpeciesAbundance(
	zone.getNucnetZone(),
	Libnucnet__Nuc__getSpeciesByName(
	  Libnucnet__Net__getNuc(
            Libnucnet__Zone__getNet( zone.getNucnetZone() )
          ),
	  Libnucnet__Reaction__Element__getName(
            reactant.getNucnetReactionElement()
          )
        )
      );

  }

  //==========================================================================
  // Iterate the products to get the reverse flow.
  //==========================================================================

  nnt::reaction_element_list_t product_list =
    nnt::make_reaction_nuclide_product_list( p_reaction );

  i_elements = product_list.size();

  d_r =
    d_reverse_rate *
    pow(
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ),
      (double) i_elements - 1.
    ) /
      Libnucnet__Reaction__getDuplicateProductFactor( p_reaction );

  BOOST_FOREACH( nnt::ReactionElement product, product_list )
  {

    d_r *=
      Libnucnet__Zone__getSpeciesAbundance(
        zone.getNucnetZone(),
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc(
            Libnucnet__Zone__getNet( zone.getNucnetZone() )
        ),
	Libnucnet__Reaction__Element__getName(
	  product.getNucnetReactionElement()
	)
      )
    );

  }

  //==========================================================================
  // Return the flows.
  //==========================================================================
  
  return std::make_pair( d_f, d_r ); 

}

//############################################################################
// compute_modified_flows_for_reaction().
//##########################################################################//

/**
 * \brief Compute the modified forward and reverse flows for a reaction.
 *
 * \param zone A Nucnet Tools zone.
 * \param p_reaction A pointer to the Libnucnet__Reaction.
 * \return A pair in which the first element is the forward flow and the
 *         second is the reverse.
 */

std::pair<double,double>
compute_modified_flows_for_reaction(
  nnt::Zone& zone,
  Libnucnet__Reaction *p_reaction
)
{

  double d_forward_rate, d_reverse_rate;
  double d_f, d_r, d_factor_forward = 0, d_factor_reverse = 0, d_mod;
  double d_factor_part;
  size_t i_elements;

  Libnucnet__Net__computeRatesForReaction(
    Libnucnet__Zone__getNet( zone.getNucnetZone() ),
    p_reaction,
    boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ),
    boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ),
    Libnucnet__Zone__getDataForUserRateFunction(
      zone.getNucnetZone(),
      Libnucnet__Reaction__getRateFunctionKey( p_reaction )
    ),
    &d_forward_rate,
    &d_reverse_rate
  );

  if( Libnucnet__Reaction__isWeak( p_reaction ) )
  {

    d_reverse_rate =
      nnt::compute_reverse_weak_rate_for_reaction(
	Libnucnet__Net__getNuc(
          Libnucnet__Zone__getNet( zone.getNucnetZone() )
        ),
	p_reaction,
	d_forward_rate,
	boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ),
	boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ),
	boost::lexical_cast<double>( zone.getProperty( nnt::s_MUEKT ) ),
	boost::lexical_cast<double>( zone.getProperty( nnt::s_MU_NUE_KT ) )
      );

  }

  //==========================================================================
  // Iterate the reactants to get the forward flow.
  //==========================================================================

  nnt::reaction_element_list_t reactant_list =
    nnt::make_reaction_nuclide_reactant_list( p_reaction );

  nnt::reaction_element_list_t mod_reactant_list =
    nnt::make_reaction_nuclide_reactant_list( p_reaction );

  i_elements = reactant_list.size();

  d_f =
    d_forward_rate *
    pow(
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ),
      (double) i_elements - 1.
    ) /
      Libnucnet__Reaction__getDuplicateReactantFactor( p_reaction );

  d_mod = d_f;

  BOOST_FOREACH( nnt::ReactionElement reactant, reactant_list )
  {

     d_f *=
       Libnucnet__Zone__getSpeciesAbundance(
	 zone.getNucnetZone(),
	 Libnucnet__Nuc__getSpeciesByName(
	   Libnucnet__Net__getNuc(
             Libnucnet__Zone__getNet( zone.getNucnetZone() )
           ),
	   Libnucnet__Reaction__Element__getName(
	     reactant.getNucnetReactionElement()
	   )
	 )
       );

     d_factor_part = d_mod;

     BOOST_FOREACH( nnt::ReactionElement mod_reactant, mod_reactant_list )
     {
 
       if(
	 mod_reactant.getNucnetReactionElement() !=
	 reactant.getNucnetReactionElement()
       )
       {

	 d_factor_part *=
	   Libnucnet__Zone__getSpeciesAbundance(
	     zone.getNucnetZone(),
	     Libnucnet__Nuc__getSpeciesByName(
	       Libnucnet__Net__getNuc(
                 Libnucnet__Zone__getNet( zone.getNucnetZone() )
               ),
	       Libnucnet__Reaction__Element__getName(
		 mod_reactant.getNucnetReactionElement()
	       )
	     )
	   );

       }

     }

     d_factor_forward +=
       d_factor_part *
         boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) );

  }

  //==========================================================================
  // Iterate the products to get the reverse flow.
  //==========================================================================

  nnt::reaction_element_list_t product_list =
    nnt::make_reaction_nuclide_product_list( p_reaction );

  nnt::reaction_element_list_t mod_product_list =
    nnt::make_reaction_nuclide_product_list( p_reaction );

  i_elements = product_list.size();

  d_r =
    d_reverse_rate *
    pow(
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ),
      (double) i_elements - 1.
    ) /
      Libnucnet__Reaction__getDuplicateProductFactor( p_reaction );

  d_mod = d_r;

  BOOST_FOREACH( nnt::ReactionElement product, product_list )
  {

     d_r *=
       Libnucnet__Zone__getSpeciesAbundance(
	 zone.getNucnetZone(),
	 Libnucnet__Nuc__getSpeciesByName(
	   Libnucnet__Net__getNuc(
             Libnucnet__Zone__getNet( zone.getNucnetZone() )
           ),
	   Libnucnet__Reaction__Element__getName(
	     product.getNucnetReactionElement()
	   )
	 )
       );

     d_factor_part = d_mod;

     BOOST_FOREACH( nnt::ReactionElement mod_product, mod_product_list )
     {
 
       if(
	 mod_product.getNucnetReactionElement() !=
	 product.getNucnetReactionElement()
       )
       {

	 d_factor_part *=
	   Libnucnet__Zone__getSpeciesAbundance(
	     zone.getNucnetZone(),
	     Libnucnet__Nuc__getSpeciesByName(
	       Libnucnet__Net__getNuc(
                 Libnucnet__Zone__getNet( zone.getNucnetZone() )
               ),
	       Libnucnet__Reaction__Element__getName(
		 mod_product.getNucnetReactionElement()
	       )
	     )
	   );

       }

     }

     d_factor_reverse +=
       d_factor_part *
         boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) );

  }

  //==========================================================================
  // Modify the flows.
  //==========================================================================
  
  if( GSL_SIGN( d_factor_reverse ) != GSL_SIGN( -d_factor_reverse ) )
  {
    d_f /= ( 1. + d_factor_forward + d_factor_reverse );

    d_r /= ( 1. + d_factor_forward + d_factor_reverse );
  }

  //==========================================================================
  // Return the flows.
  //==========================================================================
  
  return std::make_pair( d_f, d_r ); 

}

//############################################################################
// compute_forward_flow_vector().
//##########################################################################//

/**
 * \brief Compute the forward flow vector for a set of nuclei and reactions.
 *
 * \param zone A Nucnet Tools zone.
 * \param p_net_view A network view specifying the nuclei and reactions to use.
 *         the network nuclei.
 * \param s_reac_xpath A string giving an XPath expression restricting
 *         the network reactions.
 * \return A pointer to a new gsl vector containing the flows for each species.
 */

gsl_vector *
compute_forward_flow_vector(
  nnt::Zone& zone,
  Libnucnet__NetView * p_net_view
)
{

  Libnucnet__Species * p_species;
  gsl_vector *p_forward_flow;
  std::pair<double,double> flows;

  p_forward_flow =
    gsl_vector_calloc(
      Libnucnet__Nuc__getNumberOfSpecies(
	Libnucnet__Net__getNuc(
          Libnucnet__Zone__getNet( zone.getNucnetZone() )
        )
      )
    );

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list(
      Libnucnet__Net__getReac( Libnucnet__NetView__getNet( p_net_view ) )
    );
    
  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    flows =
      compute_flows_for_reaction( zone, reaction.getNucnetReaction() );

    nnt::reaction_element_list_t element_list =
      nnt::make_reaction_nuclide_reactant_list( reaction.getNucnetReaction() );

    BOOST_FOREACH( nnt::ReactionElement element, element_list )
    {

      p_species =
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc(
            Libnucnet__Zone__getNet( zone.getNucnetZone() )
          ),
          Libnucnet__Reaction__Element__getName(
            element.getNucnetReactionElement()
          )
        );

      gsl_vector_set(
        p_forward_flow,
        Libnucnet__Species__getIndex( p_species ),
        gsl_vector_get(
          p_forward_flow,
          Libnucnet__Species__getIndex( p_species )
        ) + flows.first
      );

    }

  }

  return p_forward_flow;

}

//############################################################################
// compute_reverse_flow_vector().
//##########################################################################//

/**
 * \brief Compute the reverse flow vector for a set of nuclei and reactions.
 *
 * \param zone A Nucnet Tools zone.
 * \param p_net_view A network view specifying the nuclei and reactions to use.
 * \return A pointer to a new gsl vector containing the flows for each species.
 */

gsl_vector *
compute_reverse_flow_vector(
  nnt::Zone& zone,
  Libnucnet__NetView * p_net_view
)
{

  Libnucnet__Species * p_species;
  gsl_vector * p_reverse_flow;
  std::pair<double,double> flows;

  p_reverse_flow =
    gsl_vector_calloc(
      Libnucnet__Nuc__getNumberOfSpecies(
	Libnucnet__Net__getNuc(
          Libnucnet__Zone__getNet( zone.getNucnetZone() )
        )
      )
    );

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list(
      Libnucnet__Net__getReac( Libnucnet__NetView__getNet( p_net_view ) )
    );
    
  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    flows =
      compute_flows_for_reaction( zone, reaction.getNucnetReaction() );

    nnt::reaction_element_list_t element_list =
      nnt::make_reaction_nuclide_product_list( reaction.getNucnetReaction() );

    BOOST_FOREACH( nnt::ReactionElement element, element_list )
    {

      p_species =
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc(
            Libnucnet__Zone__getNet( zone.getNucnetZone() )
          ),
          Libnucnet__Reaction__Element__getName(
            element.getNucnetReactionElement()
          )
        );

      gsl_vector_set(
        p_reverse_flow,
        Libnucnet__Species__getIndex( p_species ),
        gsl_vector_get(
          p_reverse_flow,
          Libnucnet__Species__getIndex( p_species )
        ) + flows.second
      );

    }

  }

  return p_reverse_flow;

}

//############################################################################
// compute_flow_vector().
//##########################################################################//

/**
 * \brief Compute the net flow vector for a set of nuclei and reactions.
 *
 * \param zone A Nucnet Tools zone.
 * \param p_net_view A network view specifying the nuclei and reactions to use.
 */

gsl_vector *
compute_flow_vector(
  nnt::Zone& zone,
  Libnucnet__NetView * p_net_view
)
{

  Libnucnet__Species * p_species;
  gsl_vector *p_flow;
  std::pair<double,double> flows;

  p_flow =
    gsl_vector_calloc(
      Libnucnet__Nuc__getNumberOfSpecies(
	Libnucnet__Net__getNuc(
          Libnucnet__Zone__getNet( zone.getNucnetZone() )
        )
      )
    );

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list(
      Libnucnet__Net__getReac( Libnucnet__NetView__getNet( p_net_view ) )
    );
    
  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    flows =
      compute_flows_for_reaction( zone, reaction.getNucnetReaction() );

    nnt::reaction_element_list_t reactant_list =
      nnt::make_reaction_nuclide_reactant_list(
        reaction.getNucnetReaction()
      );

    BOOST_FOREACH( nnt::ReactionElement reactant, reactant_list )
    {

      p_species =
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc(
            Libnucnet__Zone__getNet( zone.getNucnetZone() )
          ),
          Libnucnet__Reaction__Element__getName(
            reactant.getNucnetReactionElement()
          )
        );

      gsl_vector_set(
        p_flow,
        Libnucnet__Species__getIndex( p_species ),
        gsl_vector_get(
          p_flow,
          Libnucnet__Species__getIndex( p_species )
        ) + flows.first - flows.second
      );

    }

    nnt::reaction_element_list_t product_list =
      nnt::make_reaction_nuclide_product_list(
        reaction.getNucnetReaction()
      );

    BOOST_FOREACH( nnt::ReactionElement product, product_list )
    {

      p_species =
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc(
            Libnucnet__Zone__getNet( zone.getNucnetZone() )
          ),
          Libnucnet__Reaction__Element__getName(
            product.getNucnetReactionElement()
          )
        );

      gsl_vector_set(
        p_flow,
        Libnucnet__Species__getIndex( p_species ),
        gsl_vector_get(
          p_flow,
          Libnucnet__Species__getIndex( p_species )
        ) - flows.first + flows.second
      );

    }

  }

  return p_flow;

}

//############################################################################
// compute_modified_flow_vector().
//##########################################################################//

/**
 * \brief Compute the net flow vector for a set of nuclei and reactions.
 *
 * \param zone A Nucnet Tools zone.
 * \param p_net_view A network view specifying the nuclei and reactions to use.
 * \return A pointer to a new gsl vector containing the flows for each species.
 */

gsl_vector *
compute_modified_flow_vector(
  nnt::Zone& zone,
  Libnucnet__NetView * p_net_view
)
{

  Libnucnet__Species * p_species;
  gsl_vector *p_flow;
  std::pair<double,double> flows;

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list(
      Libnucnet__Net__getReac( Libnucnet__NetView__getNet( p_net_view ) )
    );
    
  p_flow =
    gsl_vector_calloc(
      Libnucnet__Nuc__getNumberOfSpecies(
	Libnucnet__Net__getNuc(
          Libnucnet__Zone__getNet( zone.getNucnetZone() )
        )
      )
    );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    flows =
      compute_modified_flows_for_reaction(
        zone,
        reaction.getNucnetReaction()
      );

    nnt::reaction_element_list_t reactant_list =
      nnt::make_reaction_nuclide_reactant_list(
        reaction.getNucnetReaction()
      );

    BOOST_FOREACH( nnt::ReactionElement reactant, reactant_list )
    {

      p_species =
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc(
            Libnucnet__Zone__getNet( zone.getNucnetZone() )
          ),
          Libnucnet__Reaction__Element__getName(
            reactant.getNucnetReactionElement()
          )
        );

      gsl_vector_set(
        p_flow,
        Libnucnet__Species__getIndex( p_species ),
        gsl_vector_get(
          p_flow,
          Libnucnet__Species__getIndex( p_species )
        ) + flows.first - flows.second
      );

    }

    nnt::reaction_element_list_t product_list =
      nnt::make_reaction_nuclide_product_list(
        reaction.getNucnetReaction()
      );

    BOOST_FOREACH( nnt::ReactionElement product, product_list )
    {

      p_species =
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc(
            Libnucnet__Zone__getNet( zone.getNucnetZone( ) )
          ),
          Libnucnet__Reaction__Element__getName(
            product.getNucnetReactionElement()
          )
        );

      gsl_vector_set(
        p_flow,
        Libnucnet__Species__getIndex( p_species ),
        gsl_vector_get(
          p_flow,
          Libnucnet__Species__getIndex( p_species )
        ) - flows.first + flows.second
      );

    }

  }

  return p_flow;

}

//############################################################################
// compute_total_flow().
//##########################################################################//

/**
 * \brief Compute the total net flow for a set of nuclei and their reactions.
 *
 * \param zone A Nucnet Tools zone.
 * \param s_type_flow A string giving the type of flow ("forward", "reverse",
 *    or "net").
 * \param A network view specifying the nuclei and reactions to use.
 * \return A double giving the sum of the net flows for the chosen reactions.
 */

double
compute_total_flow(
  nnt::Zone& zone,
  const char * s_type_flow,
  Libnucnet__NetView * p_net_view
)
{

  gsl_vector * p_flow;
  double d_result = 0;
  size_t i;

  if( strcmp( s_type_flow, nnt::s_FORWARD_FLOW ) == 0 )
    p_flow = compute_forward_flow_vector( zone, p_net_view );
  else if( strcmp( s_type_flow, nnt::s_REVERSE_FLOW ) == 0 )
    p_flow = compute_reverse_flow_vector( zone, p_net_view );
  else if( strcmp( s_type_flow, nnt::s_NET_FLOW ) == 0 )
    p_flow = compute_flow_vector( zone, p_net_view );
  else if( strcmp( s_type_flow, nnt::s_MODIFIED_NET_FLOW ) == 0 )
    p_flow = compute_modified_flow_vector( zone, p_net_view );
  else
  {
    std::cerr << "No such flow type in compute_total_flow." << std::endl;
    exit( EXIT_FAILURE );
  }

  for( i = 0; i < p_flow->size; i++ )
    d_result += gsl_vector_get( p_flow, i );

  gsl_vector_free( p_flow );

  return d_result;

}

//##############################################################################
// compute_entropy_change_rate().
//##############################################################################

/**
 * \brief Compute the entropy change rate per nucleon in a zone.
 *
 * \param zone A Nucnet Tools zone.
 * \param p_net_view A network view specifying the nuclei and reactions.
 * \return A double giving the entropy change rate per nucleon.
 */

double
compute_entropy_change_rate(
  nnt::Zone& zone,
  Libnucnet__NetView * p_net_view
)
{

  Libnucnet__Species * p_species;
  std::pair<double,double> flows;
  int i_delta_z = 0;
  double d_result = 0, d_reaction_entropy_change, d_abund;
  double d_ye, d_mu_nue_kT;
  std::vector<double> nse_corr_vector;
  coul_corr_data my_coul_corr_data;

  if(
    zone.hasProperty( nnt::s_USE_SCREENING ) &&
    zone.getProperty( nnt::s_USE_SCREENING )  == "yes"
  )
  {
  
    nnt::species_list_t species_list =
      nnt::make_species_list(
        Libnucnet__Net__getNuc( Libnucnet__NetView__getNet( p_net_view ) )
      );

    d_ye = Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 );

    my_coul_corr_data = get_coulomb_corr_data();

    BOOST_FOREACH( nnt::Species species, species_list )
    {

      nse_corr_vector.push_back( 
        my_coulomb_correction(
          species.getNucnetSpecies(),
          boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ),
          boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ),
          d_ye,
          &my_coul_corr_data
        )
      );

    }

  }

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list(
      Libnucnet__Net__getReac(
        Libnucnet__NetView__getNet( p_net_view )
      )
    );

  d_mu_nue_kT = zone.getElectronNeutrinoChemicalPotential();

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    flows =
      compute_flows_for_reaction(
        zone,
        reaction.getNucnetReaction()
      );

    i_delta_z = 0;

    if( flows.first > 1.e-100 )
    {

      d_reaction_entropy_change =
	nnt::compute_reaction_nuclear_Qvalue(
	  Libnucnet__NetView__getNet( p_net_view ),
	  reaction.getNucnetReaction(),
	  nnt::d_ELECTRON_MASS_IN_MEV
	) /
        nnt::compute_kT_in_MeV(
          boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) )
        ); 

      nnt::reaction_element_list_t reactant_list =
	nnt::make_reaction_nuclide_reactant_list(
          reaction.getNucnetReaction()
        );

      BOOST_FOREACH( nnt::ReactionElement reactant, reactant_list )
      {

	p_species =
	  Libnucnet__Nuc__getSpeciesByName(
	    Libnucnet__Net__getNuc( Libnucnet__NetView__getNet( p_net_view ) ),
	    Libnucnet__Reaction__Element__getName(
	      reactant.getNucnetReactionElement()
	    )
	  );

	i_delta_z += (int) Libnucnet__Species__getZ( p_species );

	d_abund =
	  Libnucnet__Zone__getSpeciesAbundance(
	    zone.getNucnetZone(),
	    p_species
	  );

	if( d_abund > 1.e-100 )
	{
	  d_reaction_entropy_change +=
	    log(
	      d_abund
	      /
	      Libnucnet__Species__computeQuantumAbundance(
		p_species,
		boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ),
		boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) )
	      )
	    );
	}

        if(
          zone.hasProperty( nnt::s_USE_SCREENING ) &&
          zone.getProperty( nnt::s_USE_SCREENING )  == "yes"
        )
        {
  
          d_reaction_entropy_change -=
            nse_corr_vector[ Libnucnet__Species__getIndex( p_species )]; 

        }

      }

      nnt::reaction_element_list_t product_list =
	nnt::make_reaction_nuclide_product_list(
          reaction.getNucnetReaction()
        );

      BOOST_FOREACH( nnt::ReactionElement product, product_list )
      {

	p_species =
	  Libnucnet__Nuc__getSpeciesByName(
	    Libnucnet__Net__getNuc( Libnucnet__NetView__getNet( p_net_view ) ),
	    Libnucnet__Reaction__Element__getName(
	      product.getNucnetReactionElement()
	    )
	  );

	i_delta_z -= (int) Libnucnet__Species__getZ( p_species );

	d_abund =
	  Libnucnet__Zone__getSpeciesAbundance(
	    zone.getNucnetZone(),
	    p_species
	  );

	if( d_abund > 1.e-100 )
	{
	  d_reaction_entropy_change -=
	    log(
	      d_abund
	      /
	      Libnucnet__Species__computeQuantumAbundance(
		p_species,
		boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ),
		boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) )
	      )
	    );
	}

        if(
          zone.hasProperty( nnt::s_USE_SCREENING ) &&
          zone.getProperty( nnt::s_USE_SCREENING )  == "yes"
        )
        {
  
          d_reaction_entropy_change +=
            nse_corr_vector[ Libnucnet__Species__getIndex( p_species )]; 

        }

      }

      if( i_delta_z != 0 )
      {

	d_reaction_entropy_change +=
	  (double) i_delta_z *
	  (
	     boost::lexical_cast<double>( zone.getProperty( nnt::s_MUEKT ) ) +
	     nnt::d_ELECTRON_MASS_IN_MEV /
	     nnt::compute_kT_in_MeV(
               boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) )
             )
	  );

	if( d_mu_nue_kT != GSL_NEGINF )
	{
	  d_reaction_entropy_change -= (double) i_delta_z * d_mu_nue_kT;
	}

      }

      d_result += d_reaction_entropy_change * ( flows.first - flows.second );

    }

  }

  return d_result;

}
           
//##############################################################################
// compute_energy_generation_rate_per_nucleon().
//##############################################################################

/**
 * \brief Compute the energy generation rate per nucleon in a zone.
 *
 * \param zone A Nucnet Tools zone.
 * \param p_net_view A network view specifying the nuclei and reactions.
 * \return A double giving the energy generation rate per nucleon.
 */

double
compute_energy_generation_rate_per_nucleon(
  nnt::Zone& zone,
  Libnucnet__NetView * p_net_view
)
{

  std::pair<double,double> flows;
  double d_result = 0, d_electron_energy_part;
  Libstatmech__Fermion * p_electron;
  int i_z, i_y;

  p_electron = 
    Libstatmech__Fermion__new(
      "electron",
      nnt::d_ELECTRON_MASS_IN_MEV,
      2,
      -1.
    );

  d_electron_energy_part =
    GSL_CONST_CGSM_BOLTZMANN *
    gsl_pow_2(
      boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) *
      GSL_CONST_NUM_GIGA
    ) *
    Libstatmech__Fermion__computeTemperatureDerivative(
      p_electron,
      S_CHEMICAL_POTENTIAL,
      boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) *
        GSL_CONST_NUM_GIGA,
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) *
        GSL_CONST_NUM_AVOGADRO *
        Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 ),
      NULL,
      NULL
    );

  Libstatmech__Fermion__free( p_electron );

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list(
      Libnucnet__Net__getReac(
        Libnucnet__NetView__getNet( p_net_view )
      )
    );

  for(
    nnt::reaction_list_t::iterator it = reaction_list.begin();
    it != reaction_list.end();
    it++
  )
  {

    i_z = i_y = 0;

    flows =
      compute_flows_for_reaction( zone, it->getNucnetReaction() );

    nnt::reaction_element_list_t reactant_list =
      nnt::make_reaction_nuclide_reactant_list(
        it->getNucnetReaction()
      );
        
    for(
      nnt::reaction_element_list_t::iterator iter_elem =
        reactant_list.begin();
      iter_elem != reactant_list.end();
      iter_elem++
    )
    {

      i_z -=
        (int)
        Libnucnet__Species__getZ(
          Libnucnet__Nuc__getSpeciesByName(
            Libnucnet__Net__getNuc( Libnucnet__NetView__getNet( p_net_view ) ),
            Libnucnet__Reaction__Element__getName(
              iter_elem->getNucnetReactionElement()
            )
          )
        );

      ++i_y;

    }

    nnt::reaction_element_list_t product_list =
      nnt::make_reaction_nuclide_product_list( it->getNucnetReaction() );
        
    for(
      nnt::reaction_element_list_t::iterator iter_elem =
        product_list.begin();
      iter_elem != product_list.end();
      iter_elem++
    )
    {

      i_z +=
        (int)
        Libnucnet__Species__getZ(
          Libnucnet__Nuc__getSpeciesByName(
            Libnucnet__Net__getNuc( Libnucnet__NetView__getNet( p_net_view ) ),
            Libnucnet__Reaction__Element__getName(
              iter_elem->getNucnetReactionElement()
            )
          )
        );

      --i_y;

    }

    d_result +=
      (
        (
          nnt::compute_reaction_nuclear_Qvalue(
            Libnucnet__NetView__getNet( p_net_view ),
            it->getNucnetReaction(),
            nnt::d_ELECTRON_MASS_IN_MEV
          ) *
          GSL_CONST_CGSM_ELECTRON_VOLT *
          GSL_CONST_NUM_GIGA
        )
        +
        i_y *
        (
          1.5 *
          GSL_CONST_CGSM_BOLTZMANN *
          boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) *
          GSL_CONST_NUM_GIGA
        )
        +
        i_z * 
        (
          d_electron_energy_part
          +
          nnt::d_ELECTRON_MASS_IN_MEV *
          GSL_CONST_CGSM_ELECTRON_VOLT *
          GSL_CONST_NUM_GIGA
        )
      ) * ( flows.first - flows.second );
/*
        gsl_pow_2(
          boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) *
          GSL_CONST_NUM_GIGA
        ) *
          GSL_CONST_CGSM_BOLTZMANN *
          compute_dlnG_dT(
            it->getNucnetSpecies(),
            boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) *
              GSL_CONST_NUM_GIGA
          )
*/

  }

  return d_result;

}
           
//############################################################################
// update_flow_currents().
//##########################################################################//

void
update_flow_currents( nnt::Zone& zone, nnt::Zone& flow_current_zone )
{

  update_my_rate_functions_data( zone );

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list(
      Libnucnet__Net__getReac(
        Libnucnet__NetView__getNet(
          zone.getNetView( "", "" )
        )
      )
    );

  double d_dt =
    boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    std::pair<double,double> flows =
      compute_flows_for_reaction(
        zone,
        reaction.getNucnetReaction()
      );

    double d_current = ( flows.first - flows.second ) * d_dt;

    if(
      flow_current_zone.hasProperty(
        nnt::s_FLOW_CURRENT, 
        Libnucnet__Reaction__getString( reaction.getNucnetReaction() )
      )
    )
    {
      flow_current_zone.updateProperty(
        nnt::s_FLOW_CURRENT, 
        Libnucnet__Reaction__getString( reaction.getNucnetReaction() ),
        boost::lexical_cast<std::string>(
          boost::lexical_cast<double>(
            flow_current_zone.getProperty(
              nnt::s_FLOW_CURRENT,
              Libnucnet__Reaction__getString( reaction.getNucnetReaction() )
            )
          ) + d_current
        )
      );
    }
    else
    {
      flow_current_zone.updateProperty(
        nnt::s_FLOW_CURRENT,
        Libnucnet__Reaction__getString( reaction.getNucnetReaction() ),
        boost::lexical_cast<std::string>(
          d_current
        )
      );
    }

  }

}

} // namespace user
