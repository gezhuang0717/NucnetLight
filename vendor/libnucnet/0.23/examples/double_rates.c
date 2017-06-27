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
//       Example to demonstrate how to use nucnet routines to create
//       a new Libnucnet nuclear network structure from an
//       input xml file, print out the forward and reverse rates for
//       valid reactions at the input temperature, double all reaction
//       rates and print them out again, and clear the
//       structure and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/


#include <Libnucnet.h>
#include "user_rate_functions.h"

int print_reactions( Libnucnet__Reaction *, Libnucnet__Zone * );
int double_rates( Libnucnet__Reaction *, Libnucnet__Zone * );
int
compare_reactions(
  const Libnucnet__Reaction *,
  const Libnucnet__Reaction *
);

/*##############################################################################
// main().
//############################################################################*/

int main( int argc, char * argv[] ) {

  Libnucnet__Net *p_my_network;
  Libnucnet__Zone *p_zone;
  double d_rhoe;

  /*============================================================================
  // Check input.
  //==========================================================================*/

   if ( argc < 3 || argc > 4 ) {
      fprintf(
        stderr, "\nUsage: %s filename t9 xpath\n\n", argv[0]
      );
      fprintf(
        stderr, "  filename = input network data xml filename\n\n"
      );
      fprintf(
        stderr, "  t9 = input temperature\n\n"
      );
      fprintf(
        stderr, "  xpath = reaction xpath expression (optional)\n\n"
      );

      return EXIT_FAILURE;
   }

  /*============================================================================
  // Read file and exit if not present.
  //==========================================================================*/

  p_my_network = Libnucnet__Net__new( );

  Libnucnet__Nuc__updateFromXml(
    Libnucnet__Net__getNuc( p_my_network ),
    argv[1],
    NULL
  );

  if( argv[3] )
    Libnucnet__Reac__updateFromXml(
      Libnucnet__Net__getReac( p_my_network ),
      argv[1],
      argv[3]
    );
  else
    Libnucnet__Reac__updateFromXml(
      Libnucnet__Net__getReac( p_my_network ),
      argv[1],
      NULL
    );

  /*============================================================================
  // Register user rate functions.
  //==========================================================================*/

  register_my_rate_functions( Libnucnet__Net__getReac( p_my_network ) );

  /*============================================================================
  // Create zone.
  //==========================================================================*/

  p_zone =
    Libnucnet__Zone__new(
      p_my_network,
      NULL,
      NULL,
      NULL
    );

  /*============================================================================
  // Compute rates at rho = 1 g/cc and Ye = 1.
  //==========================================================================*/

  d_rhoe = 1.;
  Libnucnet__Zone__updateDataForUserRateFunction(
    p_zone,
    CF88_WEAK_FIT,
    &d_rhoe
  );

  Libnucnet__Zone__computeRates( p_zone, atof( argv[2] ), 1. );

  /*============================================================================
  // Set the reaction compare function for alphabetically sorted
  // reactions.
  //==========================================================================*/

  Libnucnet__Reac__setReactionCompareFunction(
    Libnucnet__Net__getReac( p_my_network ),
    (Libnucnet__Reaction__compare_function) compare_reactions
  );

  /*============================================================================
  // Print rates.
  //==========================================================================*/

  printf( "\nOriginal Rates:\n" );
  printf( "\n\t\t\tReaction\t\t\t  Forward     Reverse\n" );
  printf(
    "=======================================================   =========   =========\n"
  );

  Libnucnet__Reac__iterateReactions(
    Libnucnet__Net__getReac(
      Libnucnet__Zone__getNet( p_zone )
    ),
    (Libnucnet__Reaction__iterateFunction) print_reactions,
    p_zone
  );

  printf("\n");

  /*============================================================================
  // Double the rates.  First clear the reaction compare function since
  // we don't care about the order in which we double the rates--this
  // speeds up the execution.
  //==========================================================================*/

  Libnucnet__Reac__clearReactionCompareFunction(
    Libnucnet__Net__getReac( p_my_network )
  );

  Libnucnet__Reac__iterateReactions(
    Libnucnet__Net__getReac(
      Libnucnet__Zone__getNet( p_zone )
    ),
    (Libnucnet__Reaction__iterateFunction) double_rates,
    p_zone
  );

  /*============================================================================
  // Reset the reaction compare function for alphabetically sorted
  // reactions.
  //==========================================================================*/

  Libnucnet__Reac__setReactionCompareFunction(
    Libnucnet__Net__getReac( p_my_network ),
    (Libnucnet__Reaction__compare_function) compare_reactions
  );

  /*============================================================================
  // Print new rates.
  //==========================================================================*/

  printf( "\nDoubled Rates:\n" );
  printf( "\n\t\t\tReaction\t\t\t  Forward     Reverse\n" );
  printf(
    "=======================================================   =========   =========\n"
  );

  Libnucnet__Reac__iterateReactions(
    Libnucnet__Net__getReac(
      Libnucnet__Zone__getNet( p_zone )
    ),
    (Libnucnet__Reaction__iterateFunction) print_reactions,
    p_zone
  );

  printf("\n");

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  Libnucnet__Zone__free( p_zone );
  Libnucnet__Net__free( p_my_network );

  return EXIT_SUCCESS;

}

/*##############################################################################
// print_reactions().
//############################################################################*/

int
print_reactions(
  Libnucnet__Reaction *p_reaction,
  Libnucnet__Zone *p_zone )
{

  double d_forward, d_reverse;

  if (
       Libnucnet__Net__isValidReaction(
          Libnucnet__Zone__getNet( p_zone ),
          p_reaction
       )
    ) {

      Libnucnet__Zone__getRatesForReaction(
        p_zone, p_reaction, &d_forward, &d_reverse
      );

      printf( "%-55s%12.3e%12.3e\n",
        Libnucnet__Reaction__getString( p_reaction ), d_forward, d_reverse
      );

  }

  return 1;

}

/*##############################################################################
// double_rates().
//############################################################################*/

int
double_rates( 
  Libnucnet__Reaction *p_reaction,
  Libnucnet__Zone *p_zone
)
{

  double d_forward, d_reverse;

  if (
       Libnucnet__Net__isValidReaction(
         Libnucnet__Zone__getNet( p_zone ),
         p_reaction
       )
    ) {

      Libnucnet__Zone__getRatesForReaction(
        p_zone, p_reaction, &d_forward, &d_reverse
      );

      d_forward *= 2.0;
      d_reverse *= 2.0;

      Libnucnet__Zone__updateRatesForReaction( 
        p_zone, p_reaction, d_forward, d_reverse 
      );

  }

  return 1;

}

/*##############################################################################
// compare_reactions().
//############################################################################*/

int
compare_reactions(
  const Libnucnet__Reaction *p_reaction1,
  const Libnucnet__Reaction *p_reaction2
)
{

  return
    strcmp(
      Libnucnet__Reaction__getString( p_reaction1 ),
      Libnucnet__Reaction__getString( p_reaction2 )
    );

}
