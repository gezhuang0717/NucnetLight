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
//       input xml file, print out the data for a reaction or reactions 
//       chosen by an xpath expression, and clear the structure and free
//       the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/


#include <Libnucnet__Reac.h>

void print_reaction( Libnucnet__Reaction *, void * );

int
main( int argc, char * argv[] ) {

  Libnucnet__Reac *p_my_reactions;

  /*============================================================================
  // Check input.
  //==========================================================================*/

   if ( argc != 3 ) {
      fprintf(
        stderr, "\nUsage: %s reac_filename xpath_expression\n\n", argv[0]
      );
      fprintf(
        stderr, "  reac_filename = input reactions xml filename\n\n"
      );
      fprintf(
        stderr, "  xpath_expression = reaction xpath expression\n\n"
      );
      return EXIT_FAILURE;
   }

  /*============================================================================
  // Read reaction file and exit if not present.
  //==========================================================================*/

  p_my_reactions = Libnucnet__Reac__new_from_xml( argv[1], argv[2] );

  /*============================================================================
  // Print rate data.
  //==========================================================================*/

  Libnucnet__Reac__iterateReactions(
    p_my_reactions,
    (Libnucnet__Reaction__iterateFunction) print_reaction,
    NULL
  );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  Libnucnet__Reac__free( p_my_reactions );

  return EXIT_SUCCESS;

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

  Libnucnet__Reaction__printRateData( p_reaction );

}
