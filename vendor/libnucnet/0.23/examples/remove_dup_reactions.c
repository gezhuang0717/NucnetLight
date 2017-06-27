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
//       a new Libnucnet__Reac structure from an
//       input xml file, find the duplicate reactions, remove the
//       duplicate reactions from NACRE from the reaction structure, output the
//       structure to a new xml file, and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/


#include <Libnucnet__Reac.h>

int remove_nacre_duplicate( Libnucnet__Reaction *, Libnucnet__Reac * );

int main( int argc, char * argv[] ) {

  Libnucnet__Reac *p_my_reactions, *p_my_duplicates;

  /*============================================================================
  // Check input.
  //==========================================================================*/

   if ( argc < 3 || argc > 4 ) {
      fprintf(
        stderr, "\nUsage: %s filename output_file xpath\n\n", argv[0]
      );
      fprintf(
        stderr, "  filename = input reaction xml file\n\n"
      );
      fprintf(
        stderr, "  output_file = output reaction xml file\n\n"
      );
      fprintf(
        stderr, "  xpath = reaction xpath expression (optional)\n\n"
      );

      return EXIT_FAILURE;
   }

  /*============================================================================
  // Read reactions.
  //==========================================================================*/

  if( argv[3] )
    p_my_reactions = Libnucnet__Reac__new_from_xml( argv[1], argv[3] );
  else
    p_my_reactions = Libnucnet__Reac__new_from_xml( argv[1], NULL );

  /*============================================================================
  // Print number of reactions.
  //==========================================================================*/

  fprintf(
    stdout,
    "\nThere are %lu reactions.\n\n",
    (unsigned long) Libnucnet__Reac__getNumberOfReactions( p_my_reactions )
  );

  /*============================================================================
  // Find duplicates.
  //==========================================================================*/

  p_my_duplicates = Libnucnet__Reac__getDuplicateReactions( p_my_reactions );

  /*============================================================================
  // Print number of duplicates.
  //==========================================================================*/

  fprintf(
    stdout,
    "There are %lu duplicates.\n\n",
    (unsigned long) Libnucnet__Reac__getNumberOfReactions( p_my_duplicates )
  );

  /*============================================================================
  // Remove duplicates.
  //==========================================================================*/

  Libnucnet__Reac__iterateReactions(
    p_my_duplicates,
    (Libnucnet__Reaction__iterateFunction) remove_nacre_duplicate,
    p_my_reactions
  );

  /*============================================================================
  // Print number of reactions.
  //==========================================================================*/

  fprintf(
    stdout,
    "\nThere are now %lu reactions.\n\n",
    (unsigned long) Libnucnet__Reac__getNumberOfReactions( p_my_reactions )
  );

  /*============================================================================
  // Output the cleaned up structure.
  //==========================================================================*/

  Libnucnet__Reac__writeToXmlFile( p_my_reactions, argv[2] );

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  Libnucnet__Reac__free( p_my_duplicates );
  Libnucnet__Reac__free( p_my_reactions );

  return EXIT_SUCCESS;

}

/*##############################################################################
// remove_nacre_duplicate().
//############################################################################*/

int
remove_nacre_duplicate(
  Libnucnet__Reaction *p_reaction, Libnucnet__Reac *p_reac
)
{

     printf( "%s\n", Libnucnet__Reaction__getString( p_reaction ) );

  if(
    strcmp(
      Libnucnet__Reaction__getSource( p_reaction ),
      "NACRE"
    ) == 0
  )
  {
    fprintf(
      stdout,
      "Removing %s (Data source: %s)\n",
      Libnucnet__Reaction__getString( p_reaction ),
      Libnucnet__Reaction__getSource( p_reaction )
    );
    Libnucnet__Reac__removeReaction(
      p_reac,
      Libnucnet__Reac__getReactionByString(
        p_reac, Libnucnet__Reaction__getString( p_reaction )
      )
    );
  }
  else if(
    strcmp(
      Libnucnet__Reaction__getSource(
        Libnucnet__Reaction__getParentDuplicate( p_reaction )
      ),
      "NACRE"
    ) == 0
  )
  {
    fprintf(
      stdout,
      "Removing %s (Data source: %s)\n",
      Libnucnet__Reaction__getString(
        Libnucnet__Reaction__getParentDuplicate( p_reaction )
      ),
      Libnucnet__Reaction__getSource(
        Libnucnet__Reaction__getParentDuplicate( p_reaction )
      )
    );
    Libnucnet__Reac__removeReaction(
      p_reac,
      Libnucnet__Reaction__getParentDuplicate( p_reaction )
    );
  }

  return 1;

}

