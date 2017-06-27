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
//       input xml file, print out the total number of reactions and the
//       number of reactions in a view, and clear the structure and free
//       the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libnucnet__Reac.h>

int main( int argc, char * argv[] ) {

  Libnucnet__Reac *p_my_reactions;
  Libnucnet__ReacView *p_view;

  /*============================================================================
  // Check the input.
  //==========================================================================*/

   if ( argc != 3 ) {
      fprintf(
        stderr, "\nUsage: %s reac_filename xpath_expression\n\n", argv[0]
      );
      fprintf(
        stderr, "  reac_filename = input reactions xml filename\n\n"
      );
      fprintf(
        stderr, "  xpath_expression = reaction xpath expression for view\n\n"
      );
      return EXIT_FAILURE;
   }

  /*============================================================================
  // Read reaction file
  //==========================================================================*/

  p_my_reactions = Libnucnet__Reac__new_from_xml( argv[1], NULL );

  /*============================================================================
  // Print number of reactions.
  //==========================================================================*/

  fprintf(
    stdout,
    "\nThe total number of reactions = %lu\n",
    (unsigned long) Libnucnet__Reac__getNumberOfReactions( p_my_reactions )
  );

  p_view = Libnucnet__ReacView__new( p_my_reactions, argv[2] );

  fprintf(
    stdout,
    "\nThe number of reactions in view \"%s\" = %lu\n\n",
    argv[2],
    (unsigned long)
       Libnucnet__Reac__getNumberOfReactions(
         Libnucnet__ReacView__getReac( p_view )
       )
  );

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  Libnucnet__ReacView__free( p_view );
  Libnucnet__Reac__free( p_my_reactions );

  /*============================================================================
  // Done!
  //==========================================================================*/

  return EXIT_SUCCESS;

}
