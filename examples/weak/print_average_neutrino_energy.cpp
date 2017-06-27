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
//! \brief Example code to print out average neutrino energy.
////////////////////////////////////////////////////////////////////////////////

#include <Libnucnet__Reac.h>
#include "nnt/two_d_weak_rates.h"
#include "nnt/auxiliary.h"
#include "nnt/iter.h"

#define S_XPATH_NU_E      "[product = 'neutrino_e']"
#define S_XPATH_NUBAR_E   "[product = 'anti-neutrino_e']"

int main( int argc, char * argv[] ) {

  Libnucnet__Reac *p_my_reactions; 
  Libnucnet__ReacView *p_view_nu_e, *p_view_nubar_e;
  double d_t9, d_rhoe;

  //============================================================================
  // Check input.
  //============================================================================

  if ( argc != 5 && argc!= 6 ) {
     fprintf(
       stderr, "\nUsage: %s filename t9 rho Ye xpath\n\n", argv[0]
     );
     fprintf(
       stderr, "  filename = input xml neutrino energy data filename\n\n"
     );
     fprintf(
       stderr, "  t9 = input temperature\n\n"
     );
     fprintf(
       stderr, "  rho = input density (g/cc)\n\n"
     );
     fprintf(
       stderr, "  Ye = input electron-to-baryon ratio (g/cc)\n\n"
     );
     fprintf(
       stderr, "  xpath = xpath expression (optional)\n\n"
     );

     return EXIT_FAILURE;
  }

  //============================================================================
  // Create reactions. 
  //============================================================================

  if ( argc == 5 )
    p_my_reactions = Libnucnet__Reac__new_from_xml( argv[1], NULL );
  else
    p_my_reactions = Libnucnet__Reac__new_from_xml( argv[1], argv[5] );

  //============================================================================
  // Read inputs. 
  //============================================================================

  d_t9 = atof( argv[2] );
  d_rhoe = atof( argv[3] ) * atof( argv[4] ),

  //============================================================================
  // Get reactions that emit neutrino_e.
  //============================================================================

  p_view_nu_e =
    Libnucnet__ReacView__new(
      p_my_reactions,
      S_XPATH_NU_E
    ); 

  //============================================================================
  // Get reactions that emit neutrino_e.
  //============================================================================

  p_view_nubar_e =
    Libnucnet__ReacView__new(
      p_my_reactions,
      S_XPATH_NUBAR_E
    ); 

  //============================================================================
  // Print average neutrino energy. 
  //============================================================================

  fprintf( stdout, "\n\t\tReaction\t\t\t\t\tav_e_nu(MeV)\n" );
  fprintf(
    stdout,
    "=========================================================\t==============\n"
  );

  Libnucnet__Reac__setReactionCompareFunction(
    Libnucnet__ReacView__getReac( p_view_nu_e ),
    (Libnucnet__Reaction__compare_function) nnt::compare_reactions_by_string
  );

  nnt::reaction_list_t reaction_list1 = 
    nnt::make_reaction_list( 
      Libnucnet__ReacView__getReac( p_view_nu_e )
    );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list1 )
  {

    nnt::TwoDWeakQuantity
      my_av_nu_energy(
        reaction.getNucnetReaction(),
        nnt::s_AVERAGE_ENERGY_NU_E
      );

    fprintf(
      stdout,
      "%-41s\t\t\t%12.6e\n",
      Libnucnet__Reaction__getString( reaction.getNucnetReaction() ),
      my_av_nu_energy.computeValue( d_t9, d_rhoe ).first
    );

  }

  fprintf( stdout, "\n" );

  fprintf(
    stdout,
    "Number of reactions that emit neutrino_e= %lu\n\n",
    (unsigned long)
       Libnucnet__Reac__getNumberOfReactions( 
         Libnucnet__ReacView__getReac( p_view_nu_e )
       )  
  );

  //============================================================================
  // Print average anti-neutrino energy. 
  //============================================================================

  fprintf( stdout, "\n\t\tReaction\t\t\t\t\tav_e_nubar(MeV)\n" );
  fprintf(
    stdout,
    "=========================================================\t==============\n"
  );

  Libnucnet__Reac__setReactionCompareFunction(
    Libnucnet__ReacView__getReac( p_view_nubar_e ),
    (Libnucnet__Reaction__compare_function) nnt::compare_reactions_by_string
  );

  nnt::reaction_list_t reaction_list2 = 
    nnt::make_reaction_list( 
      Libnucnet__ReacView__getReac( p_view_nubar_e )
    );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list2 )
  {

    nnt::TwoDWeakQuantity
      my_av_nu_bar_energy(
         reaction.getNucnetReaction(),
         nnt::s_AVERAGE_ENERGY_NU_E
      );

    fprintf(
      stdout,
      "%-41s\t\t\t%12.6e\n",
      Libnucnet__Reaction__getString( reaction.getNucnetReaction() ),
      my_av_nu_bar_energy.computeValue( d_t9, d_rhoe ).first
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
  // Clean up and exit.
  //============================================================================

  Libnucnet__ReacView__free( p_view_nu_e );
  Libnucnet__ReacView__free( p_view_nubar_e );
  Libnucnet__Reac__free( p_my_reactions );

  return EXIT_SUCCESS;

}

