////////////////////////////////////////////////////////////////////////////////
// This file was originally written by Bradley S. Meyer.
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//! \file
//! \brief Example code to select a valid reaction network.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <fstream>

#include <Libnucnet.h>
#include "nnt/iter.h"

#include <Libnucnet.h>

int main( int argc, char * argv[] ) {

  Libnucnet__Net *p_net, *p_new_net;
  Libnucnet__Species *p_species, *p_new_species;
  unsigned int i_z, i_a1, i_a2, i_a;
  std::ifstream my_file;

  //============================================================================
  // Check input.
  //============================================================================

   if( argc != 4 )
   {
     fprintf(
       stderr,
       "\nUsage: %s in_file mass_file out_file\n\n", argv[0]
     );
     fprintf(
       stderr, "  in_file = input xml filename\n\n"
     );
     fprintf(
       stderr, "  za_file = input network limit text filename\n\n"
     );
     fprintf(
       stderr, "  out_file = output xml filename\n\n"
     );

     exit( EXIT_FAILURE );
  }

  //============================================================================
  // Open network limit file.
  //============================================================================

  my_file.open( argv[2] );

  if( !my_file.good() )
  {
    std::cerr << "Couldn't open file " << argv[2] << "." << std::endl;
    return EXIT_FAILURE;
  }

  //============================================================================
  // Create nuclei.
  //============================================================================

  p_net =
    Libnucnet__Net__new_from_xml(
      argv[1],
      NULL,
      NULL
    );

  p_new_net = Libnucnet__Net__new( );

  //============================================================================
  // Read in nuclei.
  //============================================================================

  while( my_file >> i_z >> i_a1 >> i_a2 )
  {

    for( i_a = i_a1; i_a <= i_a2; i_a++ )
    {

      p_species =
        Libnucnet__Nuc__getSpeciesByZA(
          Libnucnet__Net__getNuc( p_net ),
          i_z,
          i_a,
          NULL
        );

      if( p_species )
      {

        p_new_species = Libnucnet__Species__copy( p_species );

        Libnucnet__Nuc__addSpecies(
          Libnucnet__Net__getNuc( p_new_net ),
          p_new_species
        );

      }

    }

  }

  my_file.close();

  //============================================================================
  // Add the valid reactions.
  //============================================================================

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list( Libnucnet__Net__getReac( p_net ) );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {
  
    if(
      Libnucnet__Net__isValidReaction(
        p_new_net,
        reaction.getNucnetReaction()
      )
    )
    {
  
      Libnucnet__Reac__addReaction(
        Libnucnet__Net__getReac( p_new_net ),
        Libnucnet__Reaction__copy( reaction.getNucnetReaction() )
      );

    }

  }

  //============================================================================
  // Write output.
  //============================================================================

  Libnucnet__Net__writeToXmlFile( p_new_net, argv[3] );

  //============================================================================
  // Clean up and done.
  //============================================================================

  Libnucnet__Net__free( p_net );
  Libnucnet__Net__free( p_new_net );
  return EXIT_SUCCESS;

}
