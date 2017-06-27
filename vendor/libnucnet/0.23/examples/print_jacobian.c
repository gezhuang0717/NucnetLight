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
//       a new Libnucnet structure from an input xml file, compute rates
//       for the input temperature and density, generate the Jacobian
//       matrix, and clear the structure and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

/*##############################################################################
// Includes.
//############################################################################*/

#include <Libnucnet.h>
#include "remove_duplicate.h"
#include "screen.h"
#include "coul_corr.h"
#include "user_rate_functions.h"

/*##############################################################################
// Define.  "no" = no screening, "yes" = screening.
//############################################################################*/

#define USE_SCREENING  "yes"

/*##############################################################################
// main().
//############################################################################*/

int
main( int argc, char * argv[] )
{

  Libnucnet *p_my_nucnet;
  Libnucnet__Zone *p_zone;
  Libnucnet__Reac *p_duplicates;
  WnMatrix *p_matrix;
  double *p_rhoe;

  /*============================================================================
  // User-supplied data structures.  See screen.h and coul_corr.h for
  // definitions.
  //==========================================================================*/

  user_screening_data my_screening_data;
  user_coul_corr_data my_corr_data;

  /*============================================================================
  // Check the input.
  //==========================================================================*/

  if ( argc!= 8 ) {
     fprintf(
       stderr,
       "\nUsage: %s in_file label_1 label_2 label_3 t9 rho out_file\n\n",
       argv[0]
     );
     fprintf(
       stderr, "  in_file = input nuclear network xml filename\n\n"
     );
     fprintf(
       stderr, "  label_1 = first label of zone\n\n"
     );
     fprintf(
       stderr, "  label_2 = second label of zone\n\n"
     );
     fprintf(
       stderr, "  label_3 = third label of zone\n\n"
     );
     fprintf(
       stderr, "  t9 = temperature (in billions of kelvins)\n\n"
     );
     fprintf(
       stderr, "  rho = density (in g/cc)\n\n"
     );
     fprintf(
       stderr, "  out_file = text file to store matrix to\n\n"
     );

     return EXIT_FAILURE;
  }

  /*============================================================================
  // Read input file.
  //==========================================================================*/

  p_my_nucnet = Libnucnet__new_from_xml( argv[1], NULL, NULL, NULL );

  /*============================================================================
  // Register user rate functions.
  //==========================================================================*/

  register_my_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

  /*============================================================================
  // Since we will allocate data for some of the user rate functions,
  // we set the rate function data deallocators.
  //==========================================================================*/

  set_user_data_deallocators(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

  /*============================================================================
  // Get zone.
  //==========================================================================*/

  p_zone =
    Libnucnet__getZoneByLabels( p_my_nucnet, argv[2], argv[3], argv[4] );

  if( !p_zone ) {
    fprintf( stderr, "No such zone!\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Remove duplicate reactions.
  //==========================================================================*/

  p_duplicates =
    Libnucnet__Reac__getDuplicateReactions(
      Libnucnet__Net__getReac(
        Libnucnet__Zone__getNet( p_zone )
      )
    );

  Libnucnet__Reac__iterateReactions(
    p_duplicates,
    (Libnucnet__Reaction__iterateFunction) remove_duplicate,
    Libnucnet__Zone__getNet( p_zone )
  );

  Libnucnet__Reac__free( p_duplicates );

  /*============================================================================
  // Set screening, if desired.
  //==========================================================================*/

  if( strcmp( USE_SCREENING, "yes" ) == 0 ) {

    /*==========================================================================
    // Set screening data and function.
    //========================================================================*/

    my_screening_data.dYe2 =
      Libnucnet__Zone__computeZMoment( p_zone, 2 );

    Libnucnet__Zone__setScreeningFunction(
      p_zone,
      (Libnucnet__Net__screening_function) my_screening_function,
      &my_screening_data
    );

    /*==========================================================================
    // Set correction function.
    //========================================================================*/

    my_corr_data.dA1 = -0.9052;
    my_corr_data.dA2 = 0.6322;

    Libnucnet__Zone__setNseCorrectionFactorFunction(
      p_zone,
      (Libnucnet__Species__nseCorrectionFactorFunction) my_coulomb_correction,
      &my_corr_data
    );

  }

  /*============================================================================
  // Set user function data.
  //==========================================================================*/

  p_rhoe = ( double * ) malloc( sizeof( double ) );

  *p_rhoe =
    atof( argv[6] ) * Libnucnet__Zone__computeZMoment( p_zone, 1 );
  
  Libnucnet__Zone__updateDataForUserRateFunction(
    p_zone,
    CF88_WEAK_FIT,
    p_rhoe
  );

  /*============================================================================
  // Compute rates.
  //==========================================================================*/

  Libnucnet__Zone__computeRates( p_zone, atof( argv[5] ), atof( argv[6] ) );

  /*============================================================================
  // Get Jacobian.
  //==========================================================================*/

  if(
      !( p_matrix = Libnucnet__Zone__computeJacobianMatrix( p_zone ) )
  ) {
      printf( "Couldn't construct Jacobian!\n" );
      return EXIT_FAILURE;
  }

  /*============================================================================
  // Write matrix to output file.
  //==========================================================================*/

  WnMatrix__writeMatrixToAsciiFile( p_matrix, argv[7], 0. );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  WnMatrix__free( p_matrix );
  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}

