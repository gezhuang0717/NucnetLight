//////////////////////////////////////////////////////////////////////////////
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
//! \brief Example code to create a latex table of two sets of reactions.
////////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <Libnucnet.h>

#include <boost/format.hpp>

#include "nnt/auxiliary.h"
#include "nnt/iter.h"

int main( int argc, char * argv[] ) {

  Libnucnet__Net * p_net;
  Libnucnet__Reac * p_reac;

  //============================================================================
  // Check input.
  //============================================================================

  if ( argc != 3 ) {
      std::cerr <<
        boost::format(
          "\nUsage: %s input_file_1 input_file_2\n\n"  
        ) % argv[0];
      std::cerr <<
        boost::format(
          "  net_file = input network data xml filename\n\n"
        );
      std::cerr <<
        boost::format(
          "  reac_file = input reaction data xml filename\n\n"
        );

      return EXIT_FAILURE;
  }

  //============================================================================
  // Read input files.
  //============================================================================

  p_net = Libnucnet__Net__new_from_xml( argv[1], NULL, NULL );

  p_reac = Libnucnet__Reac__new_from_xml( argv[2], NULL );
  
  //============================================================================
  // Beginning of table.
  //============================================================================

  std::cout << boost::format( "\\begin{table}\n" );
  std::cout << boost::format( "\\caption{}\n" );
  std::cout << boost::format( "\\centering\n" );
  std::cout << boost::format( "\\begin{tabular}{lcc}\n" );
  std::cout << boost::format( "\\hline\\hline\n" );
  std::cout << boost::format( "Reaction & Source 1 & Source 2\\\\\n" );
  std::cout << boost::format( "\\hline\n" );
  
  //============================================================================
  // Loop on reactions.
  //============================================================================

  Libnucnet__Reac__setReactionCompareFunction(
    p_reac,
    (Libnucnet__Reaction__compare_function) nnt::compare_reactions_by_string
  );

  nnt::reaction_list_t reaction_list = nnt::make_reaction_list( p_reac );

  BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
  {

    Libnucnet__Reaction * p_reaction =
      Libnucnet__Reac__getReactionByString(
        Libnucnet__Net__getReac( p_net ),
        Libnucnet__Reaction__getString( reaction.getNucnetReaction() )
    );

    if( !p_reaction )
    {
      std::cerr << "Reaction not present in first collection." << std::endl;
      return EXIT_FAILURE;
    }

    if(
      Libnucnet__Net__isValidReaction( p_net, p_reaction )
    )
    {
      std::cout <<
        boost::format( "$%s$ & %s & %s\\\\\n" ) %
        Libnucnet__Net__createValidReactionLatexString(
          p_net,
          p_reaction
        ) %
        Libnucnet__Reaction__getSource( p_reaction ) %
        Libnucnet__Reaction__getSource( reaction.getNucnetReaction() );
    }

  }

  //============================================================================
  // End of table.
  //============================================================================

  std::cout << boost::format( "\\hline\n" );
  std::cout << boost::format( "\\end{tabular}\n" );
  std::cout << boost::format( "\\label{table:my_table}\n" );
  std::cout << boost::format( "\\end{table}\n" );

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__Reac__free( p_reac );
  Libnucnet__Net__free( p_net );

  return EXIT_SUCCESS;

}

