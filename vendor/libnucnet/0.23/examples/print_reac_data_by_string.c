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
//       input xml file, print out the data for a reaction chosen by
//       its string, and clear the structure and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/


#include <Libnucnet__Reac.h>

int main( int argc, char * argv[] ) {

  Libnucnet__Reac *p_my_reactions;
  Libnucnet__Reaction *p_reaction;

  /*============================================================================
  // Check input.
  //==========================================================================*/

   if ( argc != 3 ) {
      fprintf(
        stderr, "\nUsage: %s reac_filename reac_string\n\n", argv[0]
      );
      fprintf(
        stderr, "  reac_filename = input reactions xml filename\n\n"
      );
      fprintf(
        stderr, "  reac_string = reaction string\n\n"
      );
      return EXIT_FAILURE;
   }

  /*============================================================================
  // Read reaction file and exit if not present.  For illustration, do this
  // by first creating the reaction structure and then updating the
  // data.
  //==========================================================================*/

  p_my_reactions = Libnucnet__Reac__new();

  Libnucnet__Reac__updateFromXml( p_my_reactions, argv[1], NULL );

  /*============================================================================
  // Print rate data.
  //==========================================================================*/

  p_reaction = Libnucnet__Reac__getReactionByString( p_my_reactions, argv[2] );

  if( !p_reaction ) {
    fprintf( stderr, "Reaction not found!\n" );
    return EXIT_FAILURE;
  }

  Libnucnet__Reaction__printRateData( p_reaction );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  Libnucnet__Reac__free( p_my_reactions );

  return EXIT_SUCCESS;

}

