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
//       Example to demonstrate how to use libnucnet routines to merge
//       a nuclear data xml file and a reaction data xml file into a
//       nuclear network data xml file.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/


#include <Libnucnet.h>

int main( int argc, char * argv[] ) {

  Libnucnet__Net *p_my_network;

  /*============================================================================
  // Check input.
  //==========================================================================*/

   if( argc < 4 || argc > 6 ) {
     fprintf(
       stderr,
       "\nUsage: %s nuc_file reac_file xpath_nuc xpath_reac net_file\n\n", argv[0]
     );
     fprintf(
       stderr, "  nuc_file = input nuclear data xml filename\n\n"
     );
     fprintf(
       stderr, "  reac_file = input reaction data xml filename\n\n"
     );
     fprintf(
       stderr,
       "  xpath_nuc = nuclide xpath expression (optional-required if xpath_reac present)\n\n"
     );
     fprintf(
       stderr,
       "  xpath_reac = reaction xpath expression (optional)\n\n"
     );
     fprintf(
       stderr, "  net_file = output network data xml filename\n\n"
     );

     exit( EXIT_FAILURE );
  }

  /*============================================================================
  // Create network.
  //==========================================================================*/

  p_my_network = Libnucnet__Net__new();

  /*============================================================================
  // Get input and then output to xml.
  //==========================================================================*/

  switch( argc )
  {
    case 4:
      Libnucnet__Nuc__updateFromXml(
        Libnucnet__Net__getNuc( p_my_network ), argv[1], NULL
      );
      Libnucnet__Reac__updateFromXml(
        Libnucnet__Net__getReac( p_my_network ), argv[2], NULL
      );
      Libnucnet__Net__writeToXmlFile( p_my_network, argv[3] );
      break;

    case 5:
      Libnucnet__Nuc__updateFromXml(
        Libnucnet__Net__getNuc( p_my_network ), argv[1], argv[3]
      );
      Libnucnet__Reac__updateFromXml(
        Libnucnet__Net__getReac( p_my_network ), argv[2], NULL
      );
      Libnucnet__Net__writeToXmlFile( p_my_network, argv[4] );
      break;

    case 6:
      Libnucnet__Nuc__updateFromXml(
        Libnucnet__Net__getNuc( p_my_network ), argv[1], argv[3]
      );
      Libnucnet__Reac__updateFromXml(
        Libnucnet__Net__getReac( p_my_network ), argv[2], argv[4]
      );
      Libnucnet__Net__writeToXmlFile( p_my_network, argv[5] );
      break;

    default:
      fprintf( stderr, "No such case.\n" );
      return EXIT_FAILURE;

  }

  /*============================================================================
  // Clean up and done.
  //==========================================================================*/

  Libnucnet__Net__free( p_my_network );
  return EXIT_SUCCESS;

}
