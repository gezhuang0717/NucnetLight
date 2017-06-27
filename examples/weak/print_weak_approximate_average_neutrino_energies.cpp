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
//! \brief Example code to print approximate weak rates at given t9 and rhoe.
////////////////////////////////////////////////////////////////////////////////

#include "user/aa522a25.h"

//##############################################################################
// main().
//############################################################################//

int main( int argc, char * argv[] ) {

  Libnucnet__Net *p_my_net;
  Libnucnet__ReacView *p_reac_view;
  Libstatmech__Fermion *p_electron;
  double d_eta_F;

  //============================================================================
  // Check input.
  //============================================================================

  if ( argc < 6 || argc > 8 ) {
     fprintf(
       stderr, "\nUsage: %s nuc_file t9 rho ye mu_nue_kT nuc_xpath reac_xpath\n",
       argv[0]
     );
     fprintf(
       stderr, "  nuc_file = input xml nuclear data filename\n\n"
     );
     fprintf(
       stderr, "  t9 = input T (in 10^9 K)\n\n"
     );
     fprintf(
       stderr, "  rho = input rho (in g/cc)\n\n"
     );
     fprintf(
       stderr, "  ye = input net electron to nucleon ratio\n\n"
     );
     fprintf(
       stderr, "  mu_nue_kT = input electron neutrino mu/kT\n\n"
     );
     fprintf(
       stderr, "  nuc_xpath = XPath expression for nuclei (optional--required if reac_xpath present)\n\n"
     );
     fprintf(
       stderr, "  reac_xpath = XPath expression for reactions (optional)\n\n"
     );
     return EXIT_FAILURE;
  }

  //============================================================================
  // Create network.
  //============================================================================

  p_my_net = Libnucnet__Net__new();

  if( argc == 6 )
  {
    Libnucnet__Nuc__updateFromXml(
      Libnucnet__Net__getNuc( p_my_net ),
      argv[1],
      NULL
    );
  }
  else
  {
    Libnucnet__Nuc__updateFromXml(
      Libnucnet__Net__getNuc( p_my_net ),
      argv[1],
      argv[6]
    );
  }

  //============================================================================
  // Update with approximate weak rates.
  //============================================================================

  user::aa522a25__update_net( p_my_net );

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
      atof( argv[2] ) * GSL_CONST_NUM_GIGA,
      atof( argv[3] ) * atof( argv[4] ) * GSL_CONST_NUM_AVOGADRO,
      NULL,
      NULL
    )
    +
    Libstatmech__Fermion__getRestMass( p_electron ) /
      nnt::compute_kT_in_MeV( atof( argv[2] ) );
  
  //============================================================================
  // Get reac view.
  //============================================================================

  if( argc == 8 )
    p_reac_view =
      Libnucnet__ReacView__new(
        Libnucnet__Net__getReac( p_my_net ),
        argv[7]
      );
  else
    p_reac_view =
      Libnucnet__ReacView__new(
        Libnucnet__Net__getReac( p_my_net ),
        NULL
      );
  

  //============================================================================
  // Iterate reactions.
  //============================================================================

  Libnucnet__Reac__setReactionCompareFunction(
    Libnucnet__ReacView__getReac( p_reac_view ),
    (Libnucnet__Reaction__compare_function)
       nnt::compare_reactions_by_string
  );

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list( 
      Libnucnet__ReacView__getReac( p_reac_view )
    );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {
    fprintf(
      stdout,
      "%-57s\t%g\n",
      Libnucnet__Reaction__getString( reaction.getNucnetReaction() ),
      user::aa522a25__compute_reaction_average_neutrino_energy(
	reaction.getNucnetReaction(),
        p_my_net,
        Libstatmech__Fermion__getRestMass( p_electron ),
        atof( argv[2] ),
        d_eta_F,
        atof( argv[5] )
      )
    );
  }

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libstatmech__Fermion__free( p_electron );
  Libnucnet__ReacView__free( p_reac_view );
  Libnucnet__Net__free( p_my_net );

  return EXIT_SUCCESS;

} 
