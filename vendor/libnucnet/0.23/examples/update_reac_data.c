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
//       input xml file over the web, print out those reactions and the total
//       number of reactions, update the data from a local xml file,
//       print out data about the reactions, and clear the structure and
//       free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libnucnet__Reac.h>

void print_reactions( Libnucnet__Reac * );

void print_reaction( Libnucnet__Reaction *, void * );

int main( int argc, char * argv[] ) {

  Libnucnet__Reac *p_my_reactions;

  /*============================================================================
  // Check the input.
  //==========================================================================*/

   if ( argc != 3 && argc != 4 ) {
      fprintf(
        stderr, "\nUsage: %s reac_filename new_reac_filename xpath_expression\n\n", argv[0]
      );
      fprintf(
        stderr, "  reac_filename = input reactions xml filename\n\n"
      );
      fprintf(
        stderr, "  new_reac_filename = new input reactions xml filename\n\n"
      );
      fprintf(
        stderr, "  xpath_expression = reaction xpath expression (optional)\n\n"
      );
      return EXIT_FAILURE;
   }

  /*============================================================================
  // Read reaction file
  //==========================================================================*/

  if( !argv[3] ) {
    p_my_reactions = Libnucnet__Reac__new_from_xml( argv[1], NULL );
  } else {
    p_my_reactions = Libnucnet__Reac__new_from_xml( argv[1], argv[3] );
  }

  /*============================================================================
  // Print reactions.
  //==========================================================================*/

  printf( "\nBefore reaction update:\n\n" );
  print_reactions( p_my_reactions );

  /*============================================================================
  // Update reactions.
  //==========================================================================*/

  Libnucnet__Reac__updateFromXml( p_my_reactions, argv[2], NULL );

  /*============================================================================
  // Print updated reactions.
  //==========================================================================*/

  printf( "After reaction update:\n\n" );
  print_reactions( p_my_reactions );

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  Libnucnet__Reac__free( p_my_reactions );

  /*============================================================================
  // Done!
  //==========================================================================*/

  return EXIT_SUCCESS;

}

/*##############################################################################
// print_reactions()
//############################################################################*/

void print_reactions( Libnucnet__Reac *p_my_reactions ) {

  /*============================================================================
  // Print reactions.
  //==========================================================================*/

  Libnucnet__Reac__iterateReactions(
    p_my_reactions,
    (Libnucnet__Reaction__iterateFunction)
      print_reaction,
    NULL
  );

  /*============================================================================
  // Print number of reactions.
  //==========================================================================*/

  printf(
    "\nNumber of reactions = %lu\n\n",
    (unsigned long) Libnucnet__Reac__getNumberOfReactions( p_my_reactions )
  );

}

/*##############################################################################
// print_reaction()
//############################################################################*/

void
print_reaction(
  Libnucnet__Reaction *p_reaction, void *p_user_data
)
{

  if( p_user_data ) {
    fprintf( stderr, "Routine should have no extra data" );
    exit( EXIT_FAILURE );
  }

  printf(
    "%s  (Data source: %s)\n",
    Libnucnet__Reaction__getString( p_reaction ),
    Libnucnet__Reaction__getSource( p_reaction )
  );

}
