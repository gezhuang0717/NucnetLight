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
//       input xml file, create valid network views, print out the
//       number of reactions in each view, and clear the structure and free
//       the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/


#include <Libnucnet.h>

/*##############################################################################
// Validation.  "no" = no validation, "yes" = validation.
//############################################################################*/

#define VALIDATE       "yes"

/*##############################################################################
// Defines.
//############################################################################*/

#define S_NG           "[reactant = 'n' and product = 'gamma']"
#define S_PG           "[reactant = 'h1' and product = 'gamma']"
#define S_AG           "[reactant = 'he4' and product = 'gamma']"

/*##############################################################################
// Prototypes.
//############################################################################*/

int
iterate_views(
  Libnucnet__NetView *,
  const char *,
  const char *,
  const char *,
  void *
);
int
compare_reactions(
  const Libnucnet__Reaction *,
  const Libnucnet__Reaction *
);

int
main( int argc, char * argv[] ) {

  Libnucnet__Net * p_my_network;
  Libnucnet__NetView * p_view;
  Libnucnet__Zone * p_zone;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if ( argc < 2 || argc > 4 ) {
      fprintf(
        stderr, "\nUsage: %s filename xpath_nuc xpath_reac\n\n", argv[0]
      );
      fprintf(
        stderr, "  filename = input network data xml filename\n\n"
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
  // Validate input xml file.
  //==========================================================================*/

  if( strcmp( VALIDATE, "yes" ) == 0 )
  {
    if( !Libnucnet__Net__is_valid_input_xml( argv[1] ) ) {
      fprintf( stderr, "Not a valid Libnucnet__Net input xml file!\n" );
      exit( EXIT_FAILURE );
    }
  }

  /*============================================================================
  // Read file and exit if not present.
  //==========================================================================*/

  if ( argc == 2 ) {

    if(
         !(
            p_my_network =
              Libnucnet__Net__new_from_xml( argv[1], NULL, NULL )
         )
       )
    {

      exit( EXIT_FAILURE );

    }
  
  }
  else if( argc == 3 ) {
  
    if(
         !(
            p_my_network =
              Libnucnet__Net__new_from_xml( argv[1], argv[2], NULL )
         )
       )
    {
   
      exit( EXIT_FAILURE );
  
    }

  }
  else {

    if(
         !(
            p_my_network =
              Libnucnet__Net__new_from_xml( argv[1], argv[2], argv[3] )
         )
       )
    {
   
      exit( EXIT_FAILURE );
  
    }

  }

  /*============================================================================
  // Create a zone.
  //==========================================================================*/

  p_zone = Libnucnet__Zone__new( p_my_network, "0", "0", "0" );

  /*============================================================================
  // Create a network view of all valid reactions.  Do this to save some
  // time in creating views of reaction types later.
  //==========================================================================*/

  p_view = Libnucnet__NetView__new( p_my_network, "", "" );

  /*============================================================================
  // Store views of various reaction types in zone.
  //==========================================================================*/

  Libnucnet__Zone__updateNetView(
    p_zone,
    "",
    S_NG,
    NULL,
    Libnucnet__NetView__new(
      Libnucnet__NetView__getNet( p_view ),
      "",
      S_NG
    )
  );

  Libnucnet__Zone__updateNetView(
    p_zone,
    "",
    S_PG,
    NULL,
    Libnucnet__NetView__new(
      Libnucnet__NetView__getNet( p_view ),
      "",
      S_PG
    )
  );

  Libnucnet__Zone__updateNetView(
    p_zone,
    "",
    S_AG,
    NULL,
    Libnucnet__NetView__new(
      Libnucnet__NetView__getNet( p_view ),
      "",
      S_AG
    )
  );

  /*============================================================================
  // Done with p_view, so free it.
  //==========================================================================*/

  Libnucnet__NetView__free( p_view );

  /*============================================================================
  // Set reaction compare function to output reactions in alphabetical
  // order.
  //==========================================================================*/

  Libnucnet__Reac__setReactionCompareFunction(
    Libnucnet__Net__getReac( p_my_network ),
    (Libnucnet__Reaction__compare_function) compare_reactions
  );

  /*============================================================================
  // Iterate the views.
  //==========================================================================*/

  Libnucnet__Zone__iterateNetViews(
    p_zone,
    NULL,
    NULL,
    NULL,
    (Libnucnet__NetView__iterateFunction) iterate_views,
    NULL
  );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  Libnucnet__Zone__free( p_zone );
  Libnucnet__Net__free( p_my_network );
 
  return EXIT_SUCCESS;

}

/*##############################################################################
// iterate_views().
//############################################################################*/

int
iterate_views(
  Libnucnet__NetView * p_view,
  const char * s_label1,
  const char * s_label2,
  const char * s_label3,
  void * p_data
)
{

  if( p_data || !s_label1 || s_label3 )
  {
    fprintf( stderr, "Invalid data to this routine!\n" );
    exit( EXIT_FAILURE );
  } 

  fprintf(
    stdout,
    "View %s has %lu reactions.\n",
    s_label2, 
    (unsigned long)
    Libnucnet__Reac__getNumberOfReactions(
      Libnucnet__Net__getReac( Libnucnet__NetView__getNet( p_view ) )
    )
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
