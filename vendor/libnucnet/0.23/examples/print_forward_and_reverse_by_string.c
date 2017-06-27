/*//////////////////////////////////////////////////////////////////////////////
// <file type="public">
//   <license>
//     This file was originally written by Bradley S. Meyer.
//
//     This is free software; you can redistribute it and/or modify it
//     under the terms of the GNU General Public License as published by
//     the Free Software Foundation; either version 2 of the License, or
//     (at your option) any later version.
//
//     This software is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.
//
//     Please see the src/README.txt file in this distribution for more
//     information.
//   </license>
//   <description>
//     <abstract>
//       Example to demonstrate how to use libnucnet routines to create
//       a new Libnucnet nuclear network structure from an
//       input xml file, print out the forward and reverse rates for
//       a specified reaction (chose by its string)
//       at a variety of temperatures, and clear the
//       structure and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/


#include <Libnucnet.h>
#include "user_rate_functions.h"

int main( int argc, char * argv[] ) {

  Libnucnet__Net *p_my_network;
  Libnucnet__Reaction *p_reaction;
  double d_forward, d_reverse;

  double a_t9[23] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5,
                      2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0,
                      10.0
                    };
  double d_rhoe = 1.;
  int i_temp;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if ( argc != 3 ) {
      fprintf(
        stderr, "\nUsage: %s filename string\n\n", argv[0]
      );
      fprintf(
        stderr, "  filename = input nuclear data xml filename\n\n"
      );
      fprintf(
        stderr, "  string = reaction string\n\n"
      );

      return EXIT_FAILURE;
  }

  /*============================================================================
  // Read input file.
  //==========================================================================*/

  p_my_network = Libnucnet__Net__new_from_xml( argv[1], NULL, NULL );

  /*============================================================================
  // Register user rate functions and set data.
  //==========================================================================*/

  register_my_rate_functions( Libnucnet__Net__getReac( p_my_network ) );

  /*============================================================================
  // Compute rates.
  //==========================================================================*/

  p_reaction =
    Libnucnet__Reac__getReactionByString(
      Libnucnet__Net__getReac( p_my_network ),
      argv[2]
    );

  if( !p_reaction ) {
    fprintf( stderr, "Reaction not found!\n" );
    return EXIT_FAILURE;
  }

  if( Libnucnet__Net__isValidReaction( p_my_network, p_reaction ) ) {

    printf( "\n%s\n\n", 
            Libnucnet__Reaction__getString( p_reaction )
    );

    printf( "T9(K)\t\t Rate\t\t\t Reverse\n" );
    printf( "=====\t\t============\t\t============\n" );

    for ( i_temp = 0; i_temp < 23; i_temp++ ) {

      if(
        strcmp(
          Libnucnet__Reaction__getRateFunctionKey( p_reaction ),
          CF88_WEAK_FIT
        ) == 0
      )
      {
        Libnucnet__Net__computeRatesForReaction(
          p_my_network,
          p_reaction,
          a_t9[i_temp],
          1.,
          &d_rhoe,
          &d_forward,
          &d_reverse
        );
      }
      else
      {
        Libnucnet__Net__computeRatesForReaction(
          p_my_network,
          p_reaction,
          a_t9[i_temp],
          1.,
          NULL,
          &d_forward,
          &d_reverse
        );
      }

      printf( " %1.1f\t\t%e\t\t%e\n", 
        a_t9[i_temp], d_forward, d_reverse
      );

    }

  }
  else 
  {
    fprintf( stderr, "Reaction not valid!\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  Libnucnet__Net__free( p_my_network );

  return EXIT_SUCCESS;

}
