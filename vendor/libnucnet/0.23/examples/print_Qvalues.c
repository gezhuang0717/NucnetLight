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
//       input xml file, print out the Q value for valid reactions,
//       and clear the structure and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/


#include <Libnucnet.h>

int print_valid_reaction( Libnucnet__Reaction *, Libnucnet__Net * );
int
compare_reactions(
  const Libnucnet__Reaction *, 
  const Libnucnet__Reaction *
);

int main( int argc, char * argv[] ) {

  Libnucnet__Net *p_my_network;

  /*============================================================================
  // Check input.
  //==========================================================================*/

   if( argc < 3 || argc > 5 ) {
     fprintf(
       stderr, "\nUsage: %s nuc_file reac_file xpath_nuc xpath_reac\n\n", argv[0]
     );
     fprintf(
       stderr, "  nuc_file = input nuclear data xml filename\n\n"
     );
     fprintf(
       stderr, "  reac_file = input nuclear data xml filename\n\n"
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
  // For illustration, create network from nuclear and reaction data separately
  // instead of from single net file.
  //==========================================================================*/

  p_my_network = Libnucnet__Net__new();

  if( argc > 3 ) {
    Libnucnet__Nuc__updateFromXml(
      Libnucnet__Net__getNuc( p_my_network ), argv[1], argv[3]
    );
  } else {
    Libnucnet__Nuc__updateFromXml(
      Libnucnet__Net__getNuc( p_my_network ), argv[1], NULL
    );
  }

  if( argc > 4 ) {
    Libnucnet__Reac__updateFromXml(
      Libnucnet__Net__getReac( p_my_network ), argv[2], argv[4]
    );
  } else {
    Libnucnet__Reac__updateFromXml(
      Libnucnet__Net__getReac( p_my_network ), argv[2], NULL
    );
  }

  /*============================================================================
  // Print reactions and q values.  Output in alphabetical order.
  //==========================================================================*/

  Libnucnet__Reac__setReactionCompareFunction(
    Libnucnet__Net__getReac( p_my_network ),
    (Libnucnet__Reaction__compare_function) compare_reactions
  );

  Libnucnet__Reac__iterateReactions(
    Libnucnet__Net__getReac( p_my_network ),
    (Libnucnet__Reaction__iterateFunction) print_valid_reaction,
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
  Libnucnet__Reaction *p_reaction, Libnucnet__Net *p_net 
)
{

  if( Libnucnet__Net__isValidReaction( p_net, p_reaction ) ) {

    fprintf(
      stdout,
      "Reaction:     %s\n",
      Libnucnet__Reaction__getString( p_reaction )
    );

    fprintf(
      stdout,
      "Q Value:      %f\n\n",
      Libnucnet__Net__computeReactionQValue( p_net, p_reaction )
    );

  }

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
