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
//       a new Libnucnet nuclear network structure from an
//       input xml file, check that the forward and reverse rates for
//       reactions (optionally chosen by an xpath expression)
//       are within user-specified bounds for a user-specified
//       temperature range, and clear the structure and free the
//       allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/


#include <Libnucnet.h>
#include "remove_duplicate.h"
#include "user_rate_functions.h"

#define I_N    10    /*  Number of equal logarithmic steps for temperature */

typedef struct {
  Libnucnet__Net *pNet;
  double dLowerT9;
  double dUpperT9;
  double dLowerBound;
  double dUpperBound;
} user_data;

int check_rates( Libnucnet__Reaction *, user_data * );

int main( int argc, char * argv[] ) {

  Libnucnet__Reac *p_duplicates;
  user_data *p_user_data;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if ( argc < 6 || argc > 7 ) {
      fprintf(
        stderr,
        "\nUsage: %s filename t9_lower t9_upper lower_bound upper_bound xpath\n\n",
        argv[0]
      );
      fprintf(
        stderr, "  filename = input nuclear network xml filename\n\n"
      );
      fprintf(
        stderr, "  t9_lower = lower limit on t9 range (> 0)\n\n"
      );
      fprintf(
        stderr, "  t9_upper = upper limit on t9 range\n\n"
      );
      fprintf(
        stderr, "  lower_bound = lower bound on rates\n\n"
      );
      fprintf(
        stderr, "  upper_bound = upper bound on rates\n\n"
      );
      fprintf(
        stderr, "  xpath = reaction xpath expression (optional)\n\n"
      );

      return EXIT_FAILURE;
  }

  /*============================================================================
  // Allocate memory for user structure.
  //==========================================================================*/

  p_user_data = ( user_data * ) malloc( sizeof( user_data ) );

  if( !p_user_data )
  {
    fprintf( stderr, "Couldn't allocate memory for user structure.\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Set bounds.
  //==========================================================================*/

  p_user_data->dLowerT9 = atof( argv[2] );
  p_user_data->dUpperT9 = atof( argv[3] );

  p_user_data->dLowerBound = atof( argv[4] );
  p_user_data->dUpperBound = atof( argv[5] );

  /*============================================================================
  // Check bounds.
  //==========================================================================*/

  if( p_user_data->dLowerT9 <= 0. )
  {
    fprintf( stderr, "Lower t9 must be > 0.\n" );
    return EXIT_FAILURE;
  }

  if( p_user_data->dLowerT9 >= p_user_data->dUpperT9 )
  {
    fprintf( stderr, "Lower t9 must be > upper t9.\n" );
    return EXIT_FAILURE;
  }

  if( p_user_data->dLowerBound >= p_user_data->dUpperBound )
  {
    fprintf( stderr, "Lower bound must be > upper bound.\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Read input file.
  //==========================================================================*/

  if( argc == 6 )
    p_user_data->pNet = Libnucnet__Net__new_from_xml( argv[1], NULL, NULL );
  else
    p_user_data->pNet = Libnucnet__Net__new_from_xml( argv[1], NULL, argv[6] );

  /*============================================================================
  // Register, set user rate functions, and assign data.
  //==========================================================================*/

  register_my_rate_functions(
    Libnucnet__Net__getReac( p_user_data->pNet )
  );

  /*============================================================================
  // Remove duplicate reactions.
  //==========================================================================*/

  p_duplicates =
    Libnucnet__Reac__getDuplicateReactions(
      Libnucnet__Net__getReac( p_user_data->pNet )
    );

  Libnucnet__Reac__iterateReactions(
    p_duplicates,
    (Libnucnet__Reaction__iterateFunction) remove_duplicate,
    p_user_data->pNet
  );

  Libnucnet__Reac__free( p_duplicates );

  /*============================================================================
  // Check reactions.
  //==========================================================================*/

  Libnucnet__Reac__iterateReactions(
    Libnucnet__Net__getReac( p_user_data->pNet ),
    (Libnucnet__Reaction__iterateFunction) check_rates,
    p_user_data
  );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  Libnucnet__Net__free( p_user_data->pNet );
  free( p_user_data );

  return EXIT_SUCCESS;

}

/*##############################################################################
// check_rates().
//############################################################################*/

int
check_rates( Libnucnet__Reaction *p_reaction, user_data *p_data )
{

  int i;
  double d_t9, d_forward, d_reverse, d_rhoe = 1.;

  if( Libnucnet__Net__isValidReaction( p_data->pNet, p_reaction ) )
  {

    for( i = 0; i <= I_N; i++ )
    {

      d_t9 =
        p_data->dLowerT9 *
        pow(
          ( p_data->dUpperT9 / p_data->dLowerT9 ), ( (float) i / (float) I_N )
        );

      if(
        strcmp(
          Libnucnet__Reaction__getRateFunctionKey( p_reaction ),
          CF88_WEAK_FIT
        ) == 0
      )
      {
        Libnucnet__Net__computeRatesForReaction(
          p_data->pNet,
          p_reaction,
          d_t9,
          1.,
          &d_rhoe,
          &d_forward,
          &d_reverse
        );
      }
      else
      {
        Libnucnet__Net__computeRatesForReaction(
          p_data->pNet,
          p_reaction,
          d_t9,
          1.,
          NULL,
          &d_forward,
          &d_reverse
        );
      }

      if( d_forward < p_data->dLowerBound )
        fprintf(
          stdout,
          "%s: lower bound is %g and rate at t9 = %g is %g\n", 
          Libnucnet__Reaction__getString( p_reaction ),
          p_data->dLowerBound,
          d_t9,
          d_forward
        );

      if( d_forward > p_data->dUpperBound )
        fprintf(
          stdout,
          "%s: upper bound is %g and rate at t9 = %g is %g\n", 
          Libnucnet__Reaction__getString( p_reaction ),
          p_data->dUpperBound,
          d_t9,
          d_forward
        );

      if( d_reverse < p_data->dLowerBound )
        fprintf(
          stdout,
          "%s: lower bound is %g and reverse rate at t9 = %g is %g\n", 
          Libnucnet__Reaction__getString( p_reaction ),
          p_data->dLowerBound,
          d_t9,
          d_reverse
        );

      if( d_reverse > p_data->dUpperBound )
        fprintf(
          stdout,
          "%s: upper bound is %g and reverse rate at t9 = %g is %g\n", 
          Libnucnet__Reaction__getString( p_reaction ),
          p_data->dUpperBound,
          d_t9,
          d_reverse
        );
    }

  }

  return 1;

}
