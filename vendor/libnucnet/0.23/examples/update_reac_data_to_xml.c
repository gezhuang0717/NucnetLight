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
//       input xml file, update the data from a local xml file,
//       output the updated data to xml, and clear the structure and
//       free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libnucnet__Reac.h>

int main( int argc, char * argv[] ) {

  Libnucnet__Reac *p_my_reactions;

  /*============================================================================
  // Check the input.
  //==========================================================================*/

   if ( argc < 4 || argc > 6 ) {
      fprintf(
        stderr, "\nUsage: %s reac_xml new_reac_xml updated_xml xpath_old xpath_new\n\n", argv[0]
      );
      fprintf(
        stderr, "  reac_xml = input reactions xml filename\n\n"
      );
      fprintf(
        stderr, "  new_reac_xml = new input reactions xml filename\n\n"
      );
      fprintf(
        stderr, "  update_xml = updated output reactions xml filename\n\n"
      );
      fprintf(
        stderr, "  xpath_old = reaction xpath expression on reac_xml (optional--required if xpath_new present)\n\n"
      );
      fprintf(
        stderr, "  xpath_new = reaction xpath expression on new_reac_xml (optional)\n\n"
      );
      return EXIT_FAILURE;
   }

  /*============================================================================
  // Read reaction file
  //==========================================================================*/

  if( argc < 5 ) {
    p_my_reactions = Libnucnet__Reac__new_from_xml( argv[1], NULL );
  } else {
    p_my_reactions = Libnucnet__Reac__new_from_xml( argv[1], argv[4] );
  }

  /*============================================================================
  // Update reactions.
  //==========================================================================*/

  if( argc < 6 )
    Libnucnet__Reac__updateFromXml( p_my_reactions, argv[2], NULL );
  else
    Libnucnet__Reac__updateFromXml( p_my_reactions, argv[2], argv[5] );

  /*============================================================================
  // Output to xml.
  //==========================================================================*/

  Libnucnet__Reac__writeToXmlFile( p_my_reactions, argv[3] );

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  Libnucnet__Reac__free( p_my_reactions );

  /*============================================================================
  // Done!
  //==========================================================================*/

  return EXIT_SUCCESS;

}
