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
//       a nuclear network xml file and a zone data xml file into a
//       full libnucnet data xml file.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/


#include <Libnucnet.h>

/*##############################################################################
// Validation.  "no" = no validation, "yes" = validation.
//############################################################################*/

#define VALIDATE       "yes"

int main( int argc, char * argv[] ) {

  Libnucnet *p_my_nucnet;

  /*============================================================================
  // Check input.
  //==========================================================================*/

   if( argc < 4 || argc > 7 ) {
     fprintf(
       stderr,
       "\nUsage: %s net_file zone_file xpath_nuc xpath_reac xpath_zone out_file\n\n", argv[0]
     );
     fprintf(
       stderr, "  net_file = input nuclear network xml filename\n\n"
     );
     fprintf(
       stderr, "  zone_file = input zone data xml filename\n\n"
     );
     fprintf(
       stderr,
       "  xpath_nuc = nuclide xpath expression (optional-required if xpath_reac present)\n\n"
     );
     fprintf(
       stderr,
       "  xpath_reac = reaction xpath expression (optional-required if xpath_zone present)\n\n"
     );
     fprintf(
       stderr,
       "  xpath_zone = zone xpath expression (optional)\n\n"
     );
     fprintf(
       stderr, "  out_file = output nucnet data xml filename\n\n"
     );

     exit( EXIT_FAILURE );
  }

  /*============================================================================
  // Create full libnucnet structure.
  //==========================================================================*/

  p_my_nucnet = Libnucnet__new();

  /*============================================================================
  // Validate input zone data file.
  //==========================================================================*/
  
  if( strcmp( VALIDATE, "yes" ) == 0 )
  {
    if( !Libnucnet__is_valid_zone_data_xml( argv[2] ) ) {
      fprintf( stderr, "Not valid libnucnet zone data input!\n" );
      return EXIT_FAILURE;
    }
  }
  
  /*============================================================================
  // Get input and store.
  //==========================================================================*/

  switch( argc )
  {

    case 4:
      Libnucnet__Net__updateFromXml(
        Libnucnet__getNet( p_my_nucnet ),
        argv[1],
        NULL,
        NULL
      );
      Libnucnet__assignZoneDataFromXml(
        p_my_nucnet,
        argv[2],
        NULL
      );
      Libnucnet__writeToXmlFile( p_my_nucnet, argv[3] );
      break;

    case 5:
      Libnucnet__Net__updateFromXml(
        Libnucnet__getNet( p_my_nucnet ),
        argv[1],
        argv[3],
        NULL
      );
      Libnucnet__assignZoneDataFromXml(
        p_my_nucnet,
        argv[2],
        NULL
      );
      Libnucnet__writeToXmlFile( p_my_nucnet, argv[4] );
      break;

    case 6:
      Libnucnet__Net__updateFromXml(
        Libnucnet__getNet( p_my_nucnet ),
        argv[1],
        argv[3],
        argv[4]
      );
      Libnucnet__assignZoneDataFromXml(
        p_my_nucnet,
        argv[2],
        NULL
      );
      Libnucnet__writeToXmlFile( p_my_nucnet, argv[5] );
      break;

    case 7:
      Libnucnet__Net__updateFromXml(
        Libnucnet__getNet( p_my_nucnet ),
        argv[1],
        argv[3],
        argv[4]
      );
      Libnucnet__assignZoneDataFromXml(
        p_my_nucnet,
        argv[2],
        argv[5]
      );
      Libnucnet__writeToXmlFile( p_my_nucnet, argv[6] );
      break;

    default:
      fprintf( stderr, "No such case.\n" );
      return EXIT_FAILURE;

  }

  /*============================================================================
  // Clean up and done.
  //==========================================================================*/

  Libnucnet__free( p_my_nucnet );
  return EXIT_SUCCESS;

}
