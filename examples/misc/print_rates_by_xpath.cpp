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
//! \brief Example code to compute rates at various temperatures for reactions
//          chosen by an XPath expression.
////////////////////////////////////////////////////////////////////////////////

#include <Libnucnet__Reac.h>

#include <boost/format.hpp>

#include "user/user_rate_functions.h"

#ifdef MY_USER
#include "my_user/my_rates.h"
#endif

#define D_T9_BEGIN  0.01
#define D_T9_END  10.

#define N_STEPS  100

int main( int argc, char * argv[] ) {

  Libnucnet__Reac *p_my_reactions;

  //============================================================================
  // Check input.
  //============================================================================

  if ( argc!= 3 ) {
      fprintf(
        stderr, "\nUsage: %s filename xpath\n\n", argv[0]
      );
      fprintf(
        stderr, "  filename = input nuclear data xml filename\n\n"
      );
      fprintf(
        stderr, "  xpath = xpath expression\n\n"
      );

      return EXIT_FAILURE;
  }

  //============================================================================
  // Read input file.
  //============================================================================

  p_my_reactions = Libnucnet__Reac__new_from_xml( argv[1], argv[2] );

  //============================================================================
  // Check reactions.
  //============================================================================

  if( Libnucnet__Reac__getNumberOfReactions( p_my_reactions ) == 0 )
  {
    std::cerr << "No reactions found." << std::endl;
    return EXIT_FAILURE;
  }

  //============================================================================
  // Register user-supplied rate functions and set data.
  //============================================================================

  user::register_my_rate_functions( p_my_reactions );

#ifdef MY_USER
  my_user::register_my_rate_functions( p_my_reactions );
#endif

  //============================================================================
  // Compute and print rates.
  //============================================================================

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list( p_my_reactions );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    std::cout <<
      boost::format( "\n%s\n\n" ) %
      Libnucnet__Reaction__getString( reaction.getNucnetReaction() );

    std::cout << boost::format( " T9(K)\t\t Rate\n" );
    std::cout << boost::format( "=======\t\t============\n" );

    for(
      double d_t9 = D_T9_BEGIN;
      d_t9 <= D_T9_END;
      d_t9 *=
        exp(
          ( 1. / (double) N_STEPS ) *
          log( D_T9_END / D_T9_BEGIN )
        )
    )
    {
      std::cout <<
        boost::format( "%7.4f\t\t%e\n" ) %
        d_t9 %
        Libnucnet__Reaction__computeRate(
          reaction.getNucnetReaction(),
          d_t9,
          NULL
        );
    }

  }

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__Reac__free( p_my_reactions );

  return EXIT_SUCCESS;

}

