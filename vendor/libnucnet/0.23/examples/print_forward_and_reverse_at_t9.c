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
//       input xml file, print out the forward and reverse rates at the
//       input temperature and at a density of 1 g/cc using a network view,
//       and clear the structure and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/


#include <Libnucnet.h>
#include "user_rate_functions.h"

#define I_BUFFER   256   /* Buffer for string labels */

typedef struct {
  Libnucnet__Net *pNet;
  double dT9;
} my_data;

int print_rates( Libnucnet__Reaction *, my_data * );
int
compare_reactions( const Libnucnet__Reaction *, const Libnucnet__Reaction * );

int main( int argc, char * argv[] ) {

  Libnucnet__Net *p_my_network;
  Libnucnet__NetView * p_view;
  my_data *p_my_data;
  char s_label1[I_BUFFER], s_label2[I_BUFFER];


  /*============================================================================
  // Check input.
  //==========================================================================*/

  if ( argc < 3 || argc > 5 ) {
     fprintf(
       stderr, "\nUsage: %s filename t9 xpath_nuc xpath_reac\n\n", argv[0]
     );
     fprintf(
       stderr, "  filename = input nuclear network xml filename\n\n"
     );
     fprintf(
       stderr, "  t9 = input temperature (in 10^9 K)\n\n"
     );
     fprintf(
       stderr,
       "  xpath_nuc = nuclide xpath expression (optional-required if xpath_reac present)\n\n"
     );
     fprintf(
       stderr,
       "  xpath_reac = reaction xpath expression (optional)\n\n"
     );

     exit( EXIT_FAILURE );
  }

  /*============================================================================
  // Read file and exit if not present.
  //==========================================================================*/

  if(
      !(
          p_my_network =
            Libnucnet__Net__new_from_xml( argv[1], NULL, NULL )
      )
  )
    exit( EXIT_FAILURE );

  /*============================================================================
  // Register and set user-supplied rate functions and data.
  //==========================================================================*/

  register_my_rate_functions(
    Libnucnet__Net__getReac( p_my_network )
  );

  /*============================================================================
  // Data function for reaction iteration.
  //==========================================================================*/

  p_my_data = ( my_data * ) malloc( sizeof( my_data ) );

  p_my_data->pNet = p_my_network;
  p_my_data->dT9 = atof( argv[2] );

  /*============================================================================
  // Get a network view containing all the valid reactions.
  //==========================================================================*/

  if( argc == 3 )
  {
    strcpy( s_label1, "" );
    strcpy( s_label2, "" );
  }
  else if( argc == 4 )
  {
    strcpy( s_label1, argv[3] );
    strcpy( s_label2, "" );
  }
  else
  {
    strcpy( s_label1, argv[3] );
    strcpy( s_label2, argv[4] );
  }

  p_view =
    Libnucnet__NetView__new(
      p_my_network,
      s_label1,
      s_label2
    );

  /*============================================================================
  // Print rates at input temperature and density of 1 g/cc and Ye = 1.  Ouput
  // reactions in alphabetical order.
  //==========================================================================*/

  fprintf( stdout, "\n\t\t\tReaction\t\t\t Forward     Reverse \n" );
  fprintf(
    stdout,
    "=======================================================\t ==========  ==========\n"
  );

  Libnucnet__Reac__setReactionCompareFunction(
    Libnucnet__Net__getReac(
      Libnucnet__NetView__getNet( p_view )
    ),
    (Libnucnet__Reaction__compare_function) compare_reactions
  );

  Libnucnet__Reac__iterateReactions(
    Libnucnet__Net__getReac(
      Libnucnet__NetView__getNet( p_view )
    ),
    (Libnucnet__Reaction__iterateFunction) print_rates,
    p_my_data
  );

  /*============================================================================
  // Print total number of valid reactions.
  //==========================================================================*/

  fprintf(
    stdout,
    "\n\nTotal number of valid reactions (those in network view with nuclear XPath: \"%s\" and reaction XPath: \"%s\") is %lu\n\n",
    s_label1,
    s_label2,
    (unsigned long)
      Libnucnet__Reac__getNumberOfReactions(
        Libnucnet__Net__getReac( Libnucnet__NetView__getNet( p_view ) )
    )
  );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  Libnucnet__NetView__free( p_view );
  free( p_my_data );

  Libnucnet__Net__free( p_my_network );

  return EXIT_SUCCESS;

}

/*##############################################################################
// print_rates.
//############################################################################*/

int
print_rates(
  Libnucnet__Reaction *p_reaction,
  my_data *p_my_data
)
{

  double d_forward, d_reverse, d_rhoe = 1;

  if(
    strcmp(
      Libnucnet__Reaction__getRateFunctionKey( p_reaction ),
      CF88_WEAK_FIT
    ) == 0
  )
  {
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

  fprintf(
    stdout,
    "%-55s%12.4e%12.4e\n",
    Libnucnet__Reaction__getString( p_reaction ),
    d_forward, d_reverse
  );

  return 1;

}

/*##############################################################################
// compare_reactions()
//############################################################################*/

int
compare_reactions(
  const Libnucnet__Reaction *p_reaction1,
  const Libnucnet__Reaction *p_reaction2
)
{

  return
    strcmp(
      Libnucnet__Reaction__getString( p_reaction1 ),
      Libnucnet__Reaction__getString( p_reaction2 )
    );

}
