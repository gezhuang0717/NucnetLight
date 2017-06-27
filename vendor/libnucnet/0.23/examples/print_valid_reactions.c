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
//       input xml file, print out which of the reactions in the
//       reaction data are valid (between nuclei in the network) and
//       which are invalid, and clear the structure and free
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
// Prototypes.
//############################################################################*/

int print_valid_reaction( Libnucnet__Reaction *, Libnucnet__Net * );
int print_invalid_reaction( Libnucnet__Reaction *, Libnucnet__Net * );
int
compare_reactions(
  const Libnucnet__Reaction *,
  const Libnucnet__Reaction *
);

int
main( int argc, char * argv[] ) {

  Libnucnet__Net *p_my_network;
  Libnucnet__Reac *p_reac;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if ( argc < 2 || argc > 4 ) {
      fprintf(
        stderr, "\nUsage: %s filename xpath_nuc xpath_reac\n\n", argv[0]
      );
      fprintf(
        stderr, "  filename = input nuclear data xml filename\n\n"
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
  // Print out network information.
  //==========================================================================*/

  fprintf(
    stdout,
    "\nNetwork has %lu species\n",
    (unsigned long) Libnucnet__Nuc__getNumberOfSpecies(
      Libnucnet__Net__getNuc( p_my_network )
    )
  );

  /*============================================================================
  // Get Reac structure.
  //==========================================================================*/

  p_reac = Libnucnet__Net__getReac( p_my_network );

  /*============================================================================
  // Print out reaction information.
  //==========================================================================*/

  fprintf(
    stdout,
    "\n%lu valid reactions:\n\n",
    (unsigned long) Libnucnet__Net__getNumberOfValidReactions( p_my_network )
  );

  /*----------------------------------------------------------------------------
  // Set reaction compare function to output reactions in alphabetical
  // order.
  //--------------------------------------------------------------------------*/

  Libnucnet__Reac__setReactionCompareFunction(
    p_reac,
    (Libnucnet__Reaction__compare_function) compare_reactions
  );

  /*----------------------------------------------------------------------------
  // Valid reactions.
  //--------------------------------------------------------------------------*/

  Libnucnet__Reac__iterateReactions(
    p_reac,
    (Libnucnet__Reaction__iterateFunction) print_valid_reaction,
    p_my_network
  );

  /*----------------------------------------------------------------------------
  // Invalid reactions.
  //--------------------------------------------------------------------------*/

  fprintf(
    stdout,
    "\n%lu invalid reactions:\n\n",
    (unsigned long) Libnucnet__Reac__getNumberOfReactions( p_reac )
    - (unsigned long) Libnucnet__Net__getNumberOfValidReactions( p_my_network )
  );

  Libnucnet__Reac__iterateReactions(
    p_reac,
    (Libnucnet__Reaction__iterateFunction) print_invalid_reaction,
    p_my_network
  );

  fprintf( stdout, "\n" );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  Libnucnet__Net__free( p_my_network );
 
  return EXIT_SUCCESS;

}

/*##############################################################################
// print_valid_reaction.
//############################################################################*/

int
print_valid_reaction(
  Libnucnet__Reaction *p_reaction,
  Libnucnet__Net *p_net
)
{

  if( Libnucnet__Net__isValidReaction( p_net, p_reaction ) )
    fprintf( stdout, "%s\n", Libnucnet__Reaction__getString( p_reaction ) );

  return 1;

}

/*##############################################################################
// print_invalid_reaction.
//############################################################################*/

int
print_invalid_reaction(
  Libnucnet__Reaction *p_reaction,
  Libnucnet__Net *p_net
)
{

  if( !Libnucnet__Net__isValidReaction( p_net, p_reaction ) )
    fprintf( stdout, "%s\n", Libnucnet__Reaction__getString( p_reaction ) );

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
