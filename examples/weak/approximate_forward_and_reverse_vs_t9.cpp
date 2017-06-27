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
//! \brief Example code to print weak forward and reverse rate for a
//!   reaction at a given rhoe vs. t9 based on Arcones et al., A&A 522, A25
//!   (2010).
////////////////////////////////////////////////////////////////////////////////

#include <boost/assign.hpp>

#include <Libnucnet.h>
#include <Libstatmech.h>

#include "nnt/weak_detailed_balance.h"
#include "user/aa522a25.h"

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  Libnucnet__Net * p_my_net;
  Libnucnet__Reaction * p_reaction;
  Libstatmech__Fermion * p_electron;
  std::vector<double> t9_vector =
    boost::assign::list_of(0.1)(0.2)(0.3)(0.4)(0.5)(0.6)(0.7)(0.8)(0.9)(1.0)
                          (1.5)(2.0)(2.5)(3.0)(3.5)(4.0)(4.5)(5.0)(6.0)(7.0)
                          (8.0)(9.0)(10.0);
  double d_rhoe, d_eta_F;
  double d_forward, d_reverse, d_muekT, d_mu_nue_kT;

  //============================================================================
  // Check input.
  //==========================================================================*/

  if ( argc != 6 ) {
     fprintf(
       stderr,
       "\nUsage: %s nuc_file rho Ye mu_nue_kT reaction\n\n",
       argv[0]
     );
     fprintf(
       stderr, "  nuc_file = input xml nuc data filename\n\n"
     );
     fprintf(
       stderr, "  rho = input density (g/cc)\n\n"
     );
     fprintf(
       stderr, "  Ye = input electron to baryon ratio\n\n" 
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
  // Update nuclear and reaction data.
  //==========================================================================*/

  Libnucnet__Nuc__updateFromXml(
    Libnucnet__Net__getNuc( p_my_net ), argv[1], NULL
  );

  //============================================================================
  // Update with approximate weak reaction rates. 
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
  // Loop over t9.
  //==========================================================================*/

  d_rhoe = atof( argv[2] ) * atof( argv[3] );

  d_mu_nue_kT = atof( argv[4] );

  fprintf(
    stdout,
    "\n\n%s at rhoe = %e\n\n",
    Libnucnet__Reaction__getString( p_reaction ),
    d_rhoe
  );

  fprintf( stdout, "T9(K)\t\t Forward     Reverse \n" );
  fprintf(
    stdout,
    "========== \t ==========  ==========\n"
  );



  for( size_t i_temp = 0; i_temp < t9_vector.size(); i_temp++ )
  {

    d_muekT =
      Libstatmech__Fermion__computeChemicalPotential(
        p_electron,
        t9_vector[i_temp] * GSL_CONST_NUM_GIGA,
        d_rhoe * GSL_CONST_NUM_AVOGADRO,
        NULL,
        NULL
      );

    d_eta_F =
      d_muekT
      +
      (
        Libstatmech__Fermion__getRestMass( p_electron ) /
        nnt::compute_kT_in_MeV( t9_vector[i_temp] )
      );

    d_forward =
      user::aa522a25__compute_rate(
        p_reaction,
        p_my_net,
        Libstatmech__Fermion__getRestMass( p_electron ),
        t9_vector[i_temp],
        d_eta_F,
        d_mu_nue_kT
      );

    d_reverse =
      nnt::compute_reverse_weak_rate_for_reaction(
	Libnucnet__Net__getNuc( p_my_net ),
	p_reaction,
	d_forward,
	t9_vector[i_temp],
	atof( argv[2] ),
	d_muekT,
        d_mu_nue_kT
      );

    fprintf(
      stdout,
      "%.4e\t% .4e  %.4e\n",
      t9_vector[i_temp],
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
