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
//       input xml file, print out the duplicate reactant and product
//       factors of each reaction, and clear the structure and free the
//       allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libnucnet__Reac.h>

int
compare_reactions( const Libnucnet__Reaction *, const Libnucnet__Reaction * );

void print_reaction( Libnucnet__Reaction *, void * );

int
main( int argc, char * argv[] ) {

  Libnucnet__Reac *p_my_reactions;

  /*============================================================================
  // Check input.
  //==========================================================================*/

   if ( argc != 2 && argc != 3 ) {
      fprintf(
        stderr, "\nUsage: %s reac_filename xpath_expression\n\n", argv[0]
      );
      fprintf(
        stderr, "  reac_filename = input reactions xml filename\n\n"
      );
      fprintf(
        stderr, "  xpath_expression = reaction xpath expression (optional)\n\n"
      );
      return EXIT_FAILURE;
   }

  /*============================================================================
  // Read reaction file
  //==========================================================================*/

  if( !argv[2] ) {
    p_my_reactions = Libnucnet__Reac__new_from_xml( argv[1], NULL );
  } else {
    p_my_reactions = Libnucnet__Reac__new_from_xml( argv[1], argv[2] );
  }

  /*============================================================================
  // Print reactions.
  //==========================================================================*/

  fprintf(
    stdout,
    "\nThis output shows the duplicate reactant factor (RF) and duplicate product\nfactor (PF) for each reaction in the Libnucnet__Reac structure.\n\n"
  );

  fprintf( stdout, "\n\n\t\t\tReaction\t\t\t\tRF\tPF\n\n" );

  Libnucnet__Reac__setReactionCompareFunction(
    p_my_reactions,
    (Libnucnet__Reaction__compare_function) compare_reactions
  );

  Libnucnet__Reac__iterateReactions(
    p_my_reactions,
    (Libnucnet__Reaction__iterateFunction)
       print_reaction,
    NULL
  );

  /*============================================================================
  // Do clean up.
  //==========================================================================*/

  fprintf( stdout, "\n" );

  Libnucnet__Reac__free( p_my_reactions );

  /*============================================================================
  // Done!
  //==========================================================================*/

  return EXIT_SUCCESS;

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

/*##############################################################################
// print_reaction()
//############################################################################*/

void
print_reaction( Libnucnet__Reaction *p_reaction, void *p_data )
{

  if( p_data ) {
    fprintf( stderr, "Routine should have no extra data" );
    exit( EXIT_FAILURE );
  }

  fprintf(
    stdout, 
    "%-55s  %10.1f  %6.1f\n",
    Libnucnet__Reaction__getString( p_reaction ),
    Libnucnet__Reaction__getDuplicateReactantFactor( p_reaction ),
    Libnucnet__Reaction__getDuplicateProductFactor( p_reaction )
  );

}
