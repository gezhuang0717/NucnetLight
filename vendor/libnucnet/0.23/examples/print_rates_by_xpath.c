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
//       a new Libnucnet__Reac structure of nuclear reactions from an
//       input xml file, print out the rate for a reaction or reactions
//       chosen by an xpath expression at a variety of temperatures, and clear
//       the structure and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libnucnet__Reac.h>
#include "user_rate_functions.h"

int print_rates( Libnucnet__Reaction *, void * );

int main( int argc, char * argv[] ) {

  Libnucnet__Reac *p_my_reactions;

  /*============================================================================
  // Check input.
  //==========================================================================*/

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

  /*============================================================================
  // Read input file.
  //==========================================================================*/

  p_my_reactions = Libnucnet__Reac__new_from_xml( argv[1], argv[2] );

  /*============================================================================
  // Register user-supplied rate functions and set data.
  //==========================================================================*/

  register_my_rate_functions( p_my_reactions );

  /*============================================================================
  // Compute and print rates.
  //==========================================================================*/

  Libnucnet__Reac__iterateReactions(
    p_my_reactions,
    (Libnucnet__Reaction__iterateFunction) print_rates,
    NULL
  );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  Libnucnet__Reac__free( p_my_reactions );

  return EXIT_SUCCESS;

}

/*##############################################################################
// print_rates.
//############################################################################*/

int
print_rates( Libnucnet__Reaction *p_reaction, void *p_data )
{

  double a_t9[23] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5,
                      2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0,
                      10.0
                    };
  double d_rhoe = 1;
  int i_temp;

  if( p_data ) {
    fprintf( stderr, "Routine should have no extra data" );
    exit( EXIT_FAILURE );
  }

  printf( "\n%s\n\n", 
          Libnucnet__Reaction__getString( p_reaction )
  );

  printf( "T9(K)\t\t Rate\n" );
  printf( "=====\t\t============\n" );

  for ( i_temp = 0; i_temp < 23; i_temp++ ) {

    if(
       strcmp(
         Libnucnet__Reaction__getRateFunctionKey( p_reaction ),
         CF88_WEAK_FIT
       ) == 0
    )
    {
      printf( " %1.1f\t\t%e\n", 
        a_t9[i_temp],
        Libnucnet__Reaction__computeRate(
          p_reaction,
          a_t9[i_temp],
          &d_rhoe
        )
      );
    }
    else
    {
      printf( " %1.1f\t\t%e\n", 
        a_t9[i_temp],
        Libnucnet__Reaction__computeRate(
          p_reaction,
          a_t9[i_temp],
          NULL
        )
      );
    }

  }

  return 1;

}

