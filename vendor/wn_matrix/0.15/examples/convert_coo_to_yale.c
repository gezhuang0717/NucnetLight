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
//     Please see the src/README.txt in this distribution for more copyright
//     and license information.
//   </license>
//   <description>
//     <abstract>
//       Example to demonstrate how to use wn_matrix routines to create
//       and store elements in a matrix from an input XML file, then store
//       the matrix in Yale sparse format and output the results
//       to an XML file and free the matrix and arrays.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include "WnMatrix.h"

int main( int argc, char *argv[] ) {

  WnMatrix * p_my_matrix;
  WnMatrix__Yale *p_yale;

  /*============================================================================
  // Check arguments.
  //==========================================================================*/

  if ( argc != 3 && argc != 4 ) {
    fprintf(
      stderr, "\nUsage: %s in_file yale_file \n\n", argv[0]
    );
    fprintf(
      stderr, "  in_file = input xml data file\n\n"
    );
    fprintf(
      stderr, "  yale_file = output yale xml data file\n\n"
    );
    fprintf(
      stderr,
      "  format = format code for output of matrix element values (optional--%%g if not set)\n\n"
    );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Check that input is valid.
  //==========================================================================*/

  if( !WnMatrix__is_valid_input_xml( argv[1] ) ) {
    fprintf( stderr, "Not a valid input file!\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Read in matrix from XML.
  //==========================================================================*/

  p_my_matrix = WnMatrix__new_from_xml( argv[1], NULL );

  /*============================================================================
  // Get Yale form.
  //==========================================================================*/

  p_yale = WnMatrix__getYale( p_my_matrix );

  /*============================================================================
  // Done with p_my_matrix, so free it.
  //==========================================================================*/

  WnMatrix__free( p_my_matrix );

  /*============================================================================
  // Output Yale XML file.
  //==========================================================================*/

  if( argc == 3 )
    WnMatrix__Yale__writeToXmlFile( p_yale, argv[2], NULL );
  else
    WnMatrix__Yale__writeToXmlFile( p_yale, argv[2], argv[3] );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  WnMatrix__Yale__free( p_yale );

  return EXIT_SUCCESS;

}
