//////////////////////////////////////////////////////////////////////////////
// This file was originally written by Bradley S. Meyer and Tianhong Yu.
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
//! \brief Example code to print approximate weak forward and reverse rate for a
//!   reaction at a given t9 vs. rhoe based on Arcones et al., A&A 522, A25
//!   (2010).
////////////////////////////////////////////////////////////////////////////////

#include <Libnucnet.h>
#include <Libstatmech.h>
#include "nnt/weak_detailed_balance.h"
#include "user/aa522a25.h"

//##############################################################################
// Defines.
//##############################################################################

#define I_RHO_MIN            -20       /* minimum rhoe for log10(rho/10.) */
#define I_RHO_MAX            110       /* maximum rhoe for log10(rho/10.) */

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  Libnucnet__Net * p_my_net;
  Libnucnet__Reaction * p_reaction;
  Libstatmech__Fermion * p_electron;
  double d_t9, d_rho, d_Ye, d_eta_F, d_muekT, d_mu_nue_kT;
  double d_forward, d_reverse;
  int i_rho;

  //============================================================================
  // Check input.
  //==========================================================================*/

  if ( argc != 6 ) {
     fprintf(
       stderr,
       "\nUsage: %s nuc_file t9 Ye mu_nue_kT reaction\n\n",
       argv[0]
     );
     fprintf(
       stderr, "  nuc_file = input xml nuc data filename\n\n"
     );
     fprintf(
       stderr, "  t9 = input t9 (T in 10^9 K) \n\n"
     );
     fprintf(
       stderr, "  Ye = input electron-to-nucleon ratio \n\n"
     );
     fprintf(
       stderr, "  mu_nue_kT = neutrino chemical potential / kT\n\n"
     );
     fprintf(
       stderr, "  reaction = string for reaction\n\n"
     );
     return EXIT_FAILURE;
  }

  //============================================================================
  // Create network.
  //==========================================================================*/

  p_my_net = Libnucnet__Net__new();

  //============================================================================
  // Update nuclear data.
  //==========================================================================*/

  Libnucnet__Nuc__updateFromXml(
    Libnucnet__Net__getNuc( p_my_net ), argv[1], NULL
  );

  //============================================================================
  // Update with approximate weak rates. 
  //==========================================================================*/

  user::aa522a25__update_net( p_my_net );

  //============================================================================
  // Create electron.
  //==========================================================================*/

  p_electron =
    Libstatmech__Fermion__new(
      nnt::s_ELECTRON,
      nnt::d_ELECTRON_MASS_IN_MEV,
      2,
      -1.
    );

  //============================================================================
  // Get reaction.
  //==========================================================================*/

  p_reaction =
    Libnucnet__Reac__getReactionByString(
      Libnucnet__Net__getReac( p_my_net ),
      argv[5]
    );

  if( !p_reaction )
  {
    fprintf(
      stderr,
      "Reaction %s not found.\n",
      argv[5]
    );
    return EXIT_FAILURE;
  }

  //============================================================================
  // Loop over rhoe.
  //==========================================================================*/

  d_t9 = atof( argv[2] );

  d_Ye = atof( argv[3] );

  d_mu_nue_kT = atof( argv[4] );

  fprintf(
    stdout,
    "\n\n%s at t9 = %f\n\n",
    Libnucnet__Reaction__getString( p_reaction ),
    atof( argv[2] )
  );

  fprintf( stdout, "Rhoe(g/cc)\t Forward     Reverse \n" );
  fprintf(
    stdout,
    "========== \t ==========  ==========\n"
  );

  for( i_rho = I_RHO_MIN; i_rho <= I_RHO_MAX; i_rho++ )
  {

    d_rho = pow( 10., (double) i_rho / 10. ) / d_Ye;

    d_muekT =
      Libstatmech__Fermion__computeChemicalPotential(
        p_electron,
        d_t9 * GSL_CONST_NUM_GIGA,
        d_rho * d_Ye * GSL_CONST_NUM_AVOGADRO,
        NULL,
        NULL
      );

    d_eta_F =
      d_muekT
      +
      Libstatmech__Fermion__getRestMass( p_electron ) /
      nnt::compute_kT_in_MeV( d_t9 );

    d_forward =
      user::aa522a25__compute_rate(
        p_reaction,
        p_my_net,
        Libstatmech__Fermion__getRestMass( p_electron ),
        d_t9,
        d_eta_F,
        d_mu_nue_kT
    );

    d_reverse =
      nnt::compute_reverse_weak_rate_for_reaction(
        Libnucnet__Net__getNuc( p_my_net ),
        p_reaction,
        d_forward,
        d_t9,
        d_rho,
        d_muekT,
        d_mu_nue_kT
      );

    fprintf(
      stdout,
      "%.4e\t% .4e  %.4e\n",
      d_rho * d_Ye,
      d_forward,
      d_reverse
    );

  }

  //============================================================================
  // Clean up and exit.
  //==========================================================================*/

  Libnucnet__Net__free( p_my_net );
  Libstatmech__Fermion__free( p_electron );

  return EXIT_SUCCESS;

}
