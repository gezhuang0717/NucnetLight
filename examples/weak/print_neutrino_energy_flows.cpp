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
//! \brief Example code to print out average neutrino energy and energy loss 
//!        rate.
////////////////////////////////////////////////////////////////////////////////

#include <Libnucnet__Reac.h>
#include "nnt/two_d_weak_rates.h"
#include "nnt/auxiliary.h"
#include "nnt/iter.h"

#include "user/aa522a25.h"
#include "user/weak_utilities.h"

#define S_XPATH_NU_E      "[product = 'neutrino_e']"
#define S_XPATH_NUBAR_E   "[product = 'anti-neutrino_e']"

int main( int argc, char * argv[] ) {

  Libnucnet *p_my_nucnet;
  nnt::Zone zone;
  Libnucnet__ReacView *p_view_nu_e, *p_view_nubar_e;
  Libstatmech__Fermion *p_electron;
  double d_t9, d_rhoe, d_mue_kT, d_eta_F, d_mu_nue_kT; 
  double d_energy_flow;
  double d_total_energy_loss_rate = 0.;

  //============================================================================
  // Check input.
  //============================================================================

  if ( argc != 7 && argc!= 8 ) {
     fprintf(
       stderr, "\nUsage: %s nuc_file reac_file zone_file t9 rho mu_nue_kT xpath\n\n", argv[0]
     );
     fprintf(
       stderr, "  nuc_file = input xml nuclear data filename\n\n"
     );
     fprintf(
       stderr, "  reac_file = input xml reaction data filename\n\n"
     );
     fprintf(
       stderr, "  zone_file = input xml zone data filename\n\n"
     );
     fprintf(
       stderr, "  t9 = input temperature\n\n"
     );
     fprintf(
       stderr, "  rho = input density (g/cc)\n\n"
     );
     fprintf(
       stderr, "  mu_nue_kT = electron neutrino chemical potential / kT\n\n"
     );
     fprintf(
       stderr, "  xpath = reaction xpath expression (optional)\n\n"
     );

     return EXIT_FAILURE;
  }

  //============================================================================
  // Create network and update nuclear data. 
  //============================================================================

  p_my_nucnet = Libnucnet__new();
  
  Libnucnet__Nuc__updateFromXml(
    Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
    argv[1],
    NULL
  );

  //============================================================================
  // Update reaction data. 
  //============================================================================

  if ( argc == 7 )
    Libnucnet__Reac__updateFromXml( 
      Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) ),
      argv[2], 
      NULL 
    );
  else
    Libnucnet__Reac__updateFromXml( 
      Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) ),
      argv[2], 
      argv[7] 
    );

  //============================================================================
  // Read inputs. 
  //============================================================================

  d_t9 = atof( argv[4] );
  d_mu_nue_kT = atof( argv[6] );

  //============================================================================
  // Assign zone data. 
  //============================================================================

  Libnucnet__assignZoneDataFromXml(
    p_my_nucnet,
    argv[3],
    NULL
  );

  zone.setNucnetZone( 
    Libnucnet__getZoneByLabels( p_my_nucnet, "0", "0", "0" )
  ); 

  zone.updateProperty( nnt::s_T9, argv[4] );
  zone.updateProperty( nnt::s_RHO, argv[5] );
  zone.updateProperty( nnt::s_MU_NUE_KT, argv[6] );

  d_rhoe = 
    atof( argv[5] ) *
    Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 );

  //============================================================================
  // Get reactions that emit neutrino_e and anti-neutrino_e.
  //============================================================================

  p_view_nu_e =
    Libnucnet__ReacView__new(
      Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) ), 
      S_XPATH_NU_E
    ); 

  p_view_nubar_e =
    Libnucnet__ReacView__new(
      Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) ), 
      S_XPATH_NUBAR_E
    ); 

  //============================================================================
  // Set reaction compare function. 
  //============================================================================

  Libnucnet__Reac__setReactionCompareFunction(
    Libnucnet__ReacView__getReac( p_view_nu_e ),
    (Libnucnet__Reaction__compare_function) nnt::compare_reactions_by_string
  );

  Libnucnet__Reac__setReactionCompareFunction(
    Libnucnet__ReacView__getReac( p_view_nubar_e ),
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

  d_mue_kT =
    Libstatmech__Fermion__computeChemicalPotential(
      p_electron,
      d_t9 * GSL_CONST_NUM_GIGA,
      d_rhoe * GSL_CONST_NUM_AVOGADRO,
      NULL,
      NULL
    );

  d_eta_F =
    d_mue_kT
    +
    (
      Libstatmech__Fermion__getRestMass( p_electron ) /
      nnt::compute_kT_in_MeV( d_t9 )
    );

  //============================================================================
  // Iterate reactions and print neutrino_e energy loss rate. 
  //============================================================================

  fprintf( stdout, "\n\t\t\tReaction\t\t\t\tL_nu(MeV/s)\n" );
  fprintf(
    stdout,
    "===========================================================\t==============\n"
  );

  nnt::reaction_list_t reaction_list1 = 
    nnt::make_reaction_list( 
      Libnucnet__ReacView__getReac( p_view_nu_e )
    );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list1 )
  {

    nnt::reaction_element_list_t reactant_list = 
      nnt::make_reaction_nuclide_reactant_list( reaction.getNucnetReaction() );

    nnt::reaction_element_list_t product_list = 
      nnt::make_reaction_nuclide_product_list( reaction.getNucnetReaction() );

    if( 
      reactant_list.size() == 1 &&
      product_list.size() == 1
    )
    {

      d_energy_flow = 
        user::compute_reaction_neutrino_energy_loss_rate(
          reaction.getNucnetReaction(),
          Libnucnet__Zone__getNet( zone.getNucnetZone() ),
          Libstatmech__Fermion__getRestMass( p_electron ),
          d_t9,
          d_rhoe,
          d_mu_nue_kT,
          d_eta_F
        );

      BOOST_FOREACH( nnt::ReactionElement reactant, reactant_list )
      {

        d_energy_flow *=
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

      d_total_energy_loss_rate += d_energy_flow;

    }

    fprintf(
      stdout,
      "%-57s\t%12.6e\n",
      Libnucnet__Reaction__getString( reaction.getNucnetReaction() ), 
      d_energy_flow
    );

  }

  fprintf( stdout, "\n" );

  fprintf(
    stdout,
    "Number of reactions that emit neutrino_e = %lu\n\n",
    (unsigned long)
       Libnucnet__Reac__getNumberOfReactions( 
         Libnucnet__ReacView__getReac( p_view_nu_e )
       )  
  );

  //============================================================================
  // Iterate reactions and print anti-neutrino_e energy loss rate. 
  //============================================================================

  fprintf( stdout, "\n\t\t\tReaction\t\t\t\tL_nubar(MeV/s)\n" );
  fprintf(
    stdout,
    "===========================================================\t==============\n"
  );

  nnt::reaction_list_t reaction_list2 = 
    nnt::make_reaction_list( 
      Libnucnet__ReacView__getReac( p_view_nubar_e )
    );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list2 )
  {

    nnt::reaction_element_list_t reactant_list = 
      nnt::make_reaction_nuclide_reactant_list( reaction.getNucnetReaction() );

    nnt::reaction_element_list_t product_list = 
      nnt::make_reaction_nuclide_product_list( reaction.getNucnetReaction() );

    if( 
      reactant_list.size() == 1 &&
      product_list.size() == 1
    )
    {

      d_energy_flow = 
        user::compute_reaction_neutrino_energy_loss_rate(
          reaction.getNucnetReaction(),
          Libnucnet__Zone__getNet( zone.getNucnetZone() ),
          Libstatmech__Fermion__getRestMass( p_electron ),
          d_t9,
          d_rhoe,
          d_mu_nue_kT,
          d_eta_F
        );

      BOOST_FOREACH( nnt::ReactionElement reactant, reactant_list )
      {

        d_energy_flow *=
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

      d_total_energy_loss_rate += d_energy_flow;

    }

    fprintf(
      stdout,
      "%-57s\t%12.6e\n",
      Libnucnet__Reaction__getString( reaction.getNucnetReaction() ), 
      d_energy_flow
    );

  }  

  fprintf( stdout, "\n" );

  fprintf(
    stdout,
    "Number of reactions that emit anti-neutrino_e = %lu\n\n",
    (unsigned long)
       Libnucnet__Reac__getNumberOfReactions( 
         Libnucnet__ReacView__getReac( p_view_nubar_e )
       )  
  );

  //============================================================================
  // Print the summed total neutrino energy loss rate. 
  //============================================================================

  fprintf(
    stdout,
    "The summed total neutrino_e and anti-neutrino_e energy loss rate is %g MeV/s.\n",
    d_total_energy_loss_rate
  );

  //============================================================================
  // Compute total neutrino energy loss rate. 
  //============================================================================

  fprintf(
    stdout,
    "\nThe calculated total neutrino_e and anti-neutrino_e energy loss rate is %g MeV/s.\n",
    user::compute_neutrino_energy_loss_rate( zone )
  );
      
  //============================================================================
  // Compute total neutrino entropy loss rate. 
  //============================================================================

  fprintf(
    stdout,
    "\nThe calculated total neutrino_e and anti-neutrino_e entropy loss rate is %g k_B/s.\n",
    user::compute_neutrino_entropy_loss_rate( zone )
  );
      
  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libstatmech__Fermion__free( p_electron );
  Libnucnet__ReacView__free( p_view_nu_e );
  Libnucnet__ReacView__free( p_view_nubar_e );
  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}

