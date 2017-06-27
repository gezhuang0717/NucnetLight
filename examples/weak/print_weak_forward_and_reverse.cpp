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
//! \brief Example code to print forward and reverse weak rates at
//!    a given t9, rho and Ye.
////////////////////////////////////////////////////////////////////////////////

#include <Libnucnet.h>
#include <Libstatmech.h>
#include "nnt/two_d_weak_rates.h"
#include "nnt/iter.h"
#include "nnt/weak_detailed_balance.h"
#include "user/weak_utilities.h"

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  Libnucnet__Net *p_my_net;
  Libnucnet__ReacView *p_reac_view;
  Libstatmech__Fermion * p_electron;
  double d_t9, d_rhoe, d_eta_F;
  double d_forward, d_reverse, d_muekT;

  //============================================================================
  // Check input.
  //============================================================================

  if ( argc != 7 && argc != 8 ) {
     fprintf(
       stderr, "\nUsage: %s nuc_file reac_file t9 rho Ye mu_nue_kT xpath\n\n", argv[0]
     );
     fprintf(
       stderr, "  nuc_file = input xml nuc data filename\n\n"
     );
     fprintf(
       stderr, "  reac_file = input xml reaction ft data filename\n\n"
     );
     fprintf(
       stderr, "  t9 = input temperature\n\n"
     );
     fprintf(
       stderr, "  rho = input density (g/cc)\n\n"
     );
     fprintf(
       stderr, "  Ye = input electron to baryon ratio\n\n" 
     );
     fprintf(
       stderr, "  mu_nue_kT = electron neutrino chemical potential / kT\n\n" 
     );
     fprintf(
       stderr, "  xpath = xpath expression for reactions (optional)\n\n"
     );
     return EXIT_FAILURE;
  }

  //============================================================================
  // Create network.
  //============================================================================

  p_my_net = Libnucnet__Net__new();

  //============================================================================
  // Update nuclear and reaction data.
  //============================================================================

  Libnucnet__Nuc__updateFromXml(
    Libnucnet__Net__getNuc( p_my_net ), argv[1], NULL
  );

  if( argc == 8 )
    Libnucnet__Reac__updateFromXml(
      Libnucnet__Net__getReac( p_my_net ), argv[2], argv[7]
    );
  else
    Libnucnet__Reac__updateFromXml(
      Libnucnet__Net__getReac( p_my_net ), argv[2], NULL
    );

  //============================================================================
  // Get view of two-d weak reactions.
  //============================================================================

  p_reac_view =
    Libnucnet__ReacView__new(
      Libnucnet__Net__getReac( p_my_net ),
      nnt::s_TWO_D_WEAK_XPATH
    );

  //============================================================================
  // Set reaction compare function. 
  //============================================================================

  Libnucnet__Reac__setReactionCompareFunction(
    Libnucnet__ReacView__getReac( p_reac_view ),
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

  d_rhoe = atof( argv[4] ) * atof( argv[5] );
  d_t9 = atof( argv[3] );

  d_muekT =
    Libstatmech__Fermion__computeChemicalPotential(
      p_electron,
      d_t9 * GSL_CONST_NUM_GIGA,
      d_rhoe * GSL_CONST_NUM_AVOGADRO,
      NULL,
      NULL
    );

  d_eta_F =
    d_muekT
    +
    (
      Libstatmech__Fermion__getRestMass( p_electron ) /
      nnt::compute_kT_in_MeV( d_t9 )
    );

  //============================================================================
  // Iterate reactions. 
  //============================================================================
  
  fprintf( stdout, "\n\t\t\tReaction\t\t\t Forward     Reverse \n" );
  fprintf(
    stdout,
    "=======================================================\t ==========  ==========\n");



  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list( Libnucnet__ReacView__getReac( p_reac_view ) );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    if(
      strcmp(
	Libnucnet__Reaction__getRateFunctionKey( reaction.getNucnetReaction() ),
	"two-d weak rates log10 ft"
      ) == 0
    )
    {
      d_forward =
        user::compute_weak_rate_from_log10_ft(
	  reaction.getNucnetReaction(),
          p_my_net,
          Libstatmech__Fermion__getRestMass( p_electron ),
          d_t9,
          d_rhoe,
          d_eta_F,
          atof( argv[6] )
        );
    }
    else if(
      strcmp(
	Libnucnet__Reaction__getRateFunctionKey( reaction.getNucnetReaction() ),
	"two-d weak rates"
      ) == 0
    )
    {

      d_forward =
	user::compute_two_d_weak_rate(
	  reaction.getNucnetReaction(),
	  d_t9,
          &d_rhoe
	);

    }
    else
    {
      fprintf( stderr, "Invalid weak reaction.\n" );
      return EXIT_FAILURE;
    }

  d_reverse =
    nnt::compute_reverse_weak_rate_for_reaction(
      Libnucnet__Net__getNuc( p_my_net ),
      reaction.getNucnetReaction(),
      d_forward,
      d_t9,
      d_rhoe,
      d_muekT,
      atof( argv[6] )
    );

    fprintf(
      stdout,
      "%-55s%12.4e%12.4e\n",
      Libnucnet__Reaction__getString( reaction.getNucnetReaction() ),
      d_forward,
      d_reverse
    );

  }

  //============================================================================
  // Print number of reactions.
  //============================================================================

  fprintf(
    stdout,
    "\nNumber of reactions = %lu\n\n",
    (unsigned long)
       Libnucnet__Reac__getNumberOfReactions( 
         Libnucnet__ReacView__getReac( p_reac_view )
       )
  );

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__ReacView__free( p_reac_view );
  Libnucnet__Net__free( p_my_net );
  Libstatmech__Fermion__free( p_electron );

  return EXIT_SUCCESS;

}
