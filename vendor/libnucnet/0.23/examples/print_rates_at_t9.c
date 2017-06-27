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
//       input xml file, print out the rate for each reaction at the
//       input temperature, and clear the structure and free
//       the allocated memory.  Note that the routine will interpret
//       a non-number temperature as zero.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libnucnet__Reac.h>
#include "user_rate_functions.h"

void print_rate( Libnucnet__Reaction *, double * );
int
compare_reactions(
  const Libnucnet__Reaction *,
  const Libnucnet__Reaction *
);

int main( int argc, char * argv[] ) {

  Libnucnet__Reac *p_my_reactions;
  double d_t9;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if ( argc != 3 && argc!= 4 ) {
     fprintf(
       stderr, "\nUsage: %s filename t9 xpath\n\n", argv[0]
     );
     fprintf(
       stderr, "  filename = input nuclear data xml filename\n\n"
     );
     fprintf(
       stderr, "  t9 = input temperature\n\n"
     );
     fprintf(
       stderr, "  xpath = xpath expression (optional)\n\n"
     );

     return EXIT_FAILURE;
  }

  /*============================================================================
  // Read file and exit if not present.
  //==========================================================================*/

  if ( argc == 3 ) {

    p_my_reactions = Libnucnet__Reac__new_from_xml( argv[1], NULL );

  }
  else {
  
    p_my_reactions = Libnucnet__Reac__new_from_xml( argv[1], argv[3] );

  }

  /*============================================================================
  // Register user-supplied rate functions.
  //==========================================================================*/

  register_my_rate_functions( p_my_reactions );

  /*============================================================================
  // Print rates.
  //==========================================================================*/

  fprintf( stdout, "\n\t\t\t\tReaction\t\t\tRate\n" );
  fprintf(
    stdout,
    "=========================================================\t============\n"
  );

  d_t9 = atof( argv[2] );

  Libnucnet__Reac__setReactionCompareFunction(
    p_my_reactions,
    (Libnucnet__Reaction__compare_function) compare_reactions
  );

  Libnucnet__Reac__iterateReactions(
    p_my_reactions,
    (Libnucnet__Reaction__iterateFunction) print_rate,
    &d_t9
  );

  fprintf( stdout, "\n" );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  Libnucnet__Reac__free( p_my_reactions );

  return EXIT_SUCCESS;

}

/*##############################################################################
// print_rate.
//############################################################################*/

void
print_rate( Libnucnet__Reaction *p_reaction, double *p_t9 )
{

  double d_rhoe = 1;

  if(
    strcmp(
      Libnucnet__Reaction__getRateFunctionKey( p_reaction ),
      CF88_WEAK_FIT
    ) == 0
  )
  {
    fprintf(
      stdout,
      "%-57s\t%12.6e\n",
      Libnucnet__Reaction__getString( p_reaction ),
      Libnucnet__Reaction__computeRate(
        p_reaction,
        *p_t9,
        &d_rhoe
      )
    );
  }
  else
  {
    fprintf(
      stdout,
      "%-57s\t%12.6e\n",
      Libnucnet__Reaction__getString( p_reaction ),
      Libnucnet__Reaction__computeRate(
        p_reaction,
        *p_t9,
        NULL
      )
    );
  }

}


/*##############################################################################
// compare_reactions()
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
