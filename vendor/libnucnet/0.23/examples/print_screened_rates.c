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
//       a new Libnucnet structure from an input xml file, 
//       supply a user-defined screening function and data, compute
//       rates for a set of reactions chosen by an xpath expression,
//       print out the reaction rates and the screening factor,
//       and clear the structure and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libnucnet.h>
#include "screen.h"
#include "coul_corr.h"
#include "user_rate_functions.h"

/*##############################################################################
// Structure for passing data to print_reaction.
//############################################################################*/

typedef struct
{
  Libnucnet__Net *pNet;
  double dT9, dRho, dYe;
  Libnucnet__Net__screening_function pfScreeningFunction;
  user_screening_data *pScreeningData;
  Libnucnet__Species__nseCorrectionFactorFunction
    pfNseCorrectionFactorFunction;
  user_coul_corr_data *pCoulCorrData;
} extra_data;

/*##############################################################################
// Prototype.
//############################################################################*/

int
print_reaction( Libnucnet__Reaction *, extra_data * );

/*##############################################################################
// main()
//############################################################################*/

int main( int argc, char * argv[] ) {

  extra_data my_data;

  /*============================================================================
  // User-supplied data.  Type defined in screen.h.
  //==========================================================================*/

  user_screening_data my_screening_data;
  user_coul_corr_data my_corr_data;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if ( argc!= 7 ) {
      fprintf(
        stderr, "\nUsage: %s net_file t9 rho Ye Ye2 xpath_reaction\n\n", argv[0]
      );
      fprintf(
        stderr, "  net_file = input nuclear network xml filename\n\n"
      );
      fprintf(
        stderr, "  t9 = input temperature (in 10^9 K)\n\n"
      );
      fprintf(
        stderr, "  rho = input mass density (in g/cc)\n\n"
      );
      fprintf(
        stderr, "  Ye = input electron-to-baryon ratio\n\n"
      );
      fprintf(
        stderr, "  Ye2 = input second abundance moment wrt Z\n\n"
      );
      fprintf(
        stderr, "  xpath_reaction = reaction xpath expression\n\n"
      );

      return EXIT_FAILURE;
  }

  /*============================================================================
  // Read input file.
  //==========================================================================*/

  my_data.pNet =
    Libnucnet__Net__new_from_xml(
      argv[1], NULL, argv[6]
    );

  /*============================================================================
  // Register user-supplied functions.
  //==========================================================================*/

  register_my_rate_functions( Libnucnet__Net__getReac( my_data.pNet ) );

  /*============================================================================
  // Set screening data.  The default data to the screening function are
  // t9, rho, and Ye.  To pass in extra data, use your own data structure.
  // Here we pass in Ye2 = sum of Z^2 * Y_i, that is, the second moment of
  // the abundances with respect to Z.  We use our own data structure
  // named my_screening_data and of the type we defined as user_screening_data.
  //==========================================================================*/

  my_screening_data.dYe2 = atof( argv[5] );

  /*============================================================================
  // Set NSE correction data.  The default data to the correction function are
  // t9, rho, and Ye.  To pass in extra data, use your own data structure.
  // Although the parameters A1 and A2 could be set in the correction function
  // itself, we pass them in here through our own data structure for
  // illustration.
  //==========================================================================*/

  my_corr_data.dA1 = -0.9052;
  my_corr_data.dA2 =  0.6322;

  /*============================================================================
  // Set extra data.
  //==========================================================================*/

  my_data.dT9 = atof( argv[2] );
  my_data.dRho = atof( argv[3] );
  my_data.dYe = atof( argv[4] );
  my_data.pfScreeningFunction =
    (Libnucnet__Net__screening_function) my_screening_function;
  my_data.pScreeningData = &my_screening_data;
  my_data.pfNseCorrectionFactorFunction =
    (Libnucnet__Species__nseCorrectionFactorFunction) my_coulomb_correction;
  my_data.pCoulCorrData = &my_corr_data;

  /*============================================================================
  // Loop over the selected (and valid) reactions.
  //==========================================================================*/

  Libnucnet__Reac__iterateReactions(
    Libnucnet__Net__getReac( my_data.pNet ),
    (Libnucnet__Reaction__iterateFunction) print_reaction,
    &my_data
  );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  Libnucnet__Net__free( my_data.pNet );
  return EXIT_SUCCESS;

}

/*##############################################################################
// print_reaction()
//############################################################################*/

int
print_reaction( Libnucnet__Reaction *p_reaction, extra_data *p_my_data )
{

  double d_forward, d_reverse, d_screen, d_correction, d_rhoe;

  if( Libnucnet__Net__isValidReaction( p_my_data->pNet, p_reaction ) ) {

    if(
      strcmp(
        Libnucnet__Reaction__getRateFunctionKey( p_reaction ),
        CF88_WEAK_FIT
      ) == 0
    )
    {
      d_rhoe = p_my_data->dRho * p_my_data->dYe;
      Libnucnet__Net__computeRatesForReaction(
        p_my_data->pNet,
        p_reaction,
        p_my_data->dT9,
        1.,
        &d_rhoe,
        &d_forward,
        &d_reverse
      );
    }
    else
    {
      Libnucnet__Net__computeRatesForReaction(
        p_my_data->pNet,
        p_reaction,
        p_my_data->dT9,
        1.,
        NULL,
        &d_forward,
        &d_reverse
       );
    }

    d_screen =
      Libnucnet__Net__computeScreeningFactorForReaction(
        p_my_data->pNet,
        p_reaction,
        p_my_data->dT9,
        p_my_data->dRho,
        p_my_data->dYe,
        p_my_data->pfScreeningFunction,
        p_my_data->pScreeningData
      );

    d_correction =
      Libnucnet__Net__computeReverseRatioCorrectionFactorForReaction(
        p_my_data->pNet,
        p_reaction,
        p_my_data->dT9,
        p_my_data->dRho,
        p_my_data->dYe,
        (Libnucnet__Species__nseCorrectionFactorFunction)
           p_my_data->pfNseCorrectionFactorFunction,
        p_my_data->pCoulCorrData
      );

    fprintf( stdout, "\n%s:\n", Libnucnet__Reaction__getString( p_reaction ) );
    fprintf( stdout, "  forward = %e\n", d_forward );
    fprintf( stdout, "  reverse = %e\n", d_reverse );
    fprintf( stdout, "  forward screening factor = %e\n", d_screen );
    fprintf(
      stdout,
      "  reverse ratio correction factor = %e\n",
      d_correction
    );
    fprintf( stdout, "  new forward rate = %e\n", d_screen * d_forward );
    fprintf(
      stdout,
      "  new reverse rate = %e\n",
      d_screen * d_correction * d_reverse
    );

  }

  return 1;

}
