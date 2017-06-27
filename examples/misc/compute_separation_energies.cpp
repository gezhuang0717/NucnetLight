////////////////////////////////////////////////////////////////////////////////
// This file was originally written by Tianhong Yu and Bradley S. Meyer.
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
//! \brief Example code to print out separation energies.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <iostream>

#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

#include <nnt/auxiliary.h>
#include <Libnucnet.h>

//##############################################################################
// Defines.
//##############################################################################

#define B_OUTPUT_ALL_SPECIES   false     // Change to true to print all species 
                                         // with S_BLANK for species that don't
                                         // have all the necessary data to
                                         // compute the separation energy.

#define S_BLANK               "-------"  // String to print for species for
                                         // which separation energy can't be
                                         // computed.  Replace  ------- with
                                         // desired string, such as 9999.99

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  Libnucnet__Nuc * p_nuc;
  Libnucnet__Species * p_sep, * p_species2;
  int delta;

  //============================================================================
  // Check input.
  //============================================================================

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {
    std::cout << std::endl;
    std::cout << argv[0] <<
      " ../../data_pub/my_net.xml h1 2 \"[a = 1 or a - z = 45]\"" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc < 4 || argc > 5 )
  {
    fprintf(
      stderr,
      "\nUsage: %s xml species delta\n\n",
      argv[0]
    );
    fprintf(
      stderr, "  xml = input xml filename with nuclear data\n\n"
    );
    fprintf(
      stderr,
      "  species = species for separation energy\n\n"
    );
    fprintf(
      stderr, "  delta = number of nucleons for separation\n\n"
    );
    fprintf(
      stderr, "  nuc_xpath = XPath to select nuclides (optional)\n\n"
    );
    std::cout << "For an example usage, type " << std::endl << std::endl;
    std::cout << argv[0] << " --example" << std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  //============================================================================
  // Read and store input.
  //============================================================================

  if( argc == 4 )
    p_nuc = Libnucnet__Nuc__new_from_xml( argv[1], NULL );
  else
    p_nuc = Libnucnet__Nuc__new_from_xml( argv[1], argv[4] );

  //============================================================================
  // Get nucleon.
  //============================================================================

  p_sep = Libnucnet__Nuc__getSpeciesByName( p_nuc, argv[2] );

  if( !p_sep )
  {
    std::cerr << "Species " << argv[2] <<
                 " not present in collection.  Note: use \"h1\" for p and \"h2\" for d." << std::endl;
    return EXIT_FAILURE;
  }

  //============================================================================
  // Get delta.
  //============================================================================

  delta = boost::lexical_cast<int>( argv[3] );

  if( delta < 0 )
  {
    std::cerr << "Delta must be >= 0." << std::endl;
    return EXIT_FAILURE;
  }

  //============================================================================
  // Print header.
  //============================================================================

  std::cout <<
    boost::format( "    Z      N      A      S_%s_%s (MeV)\n" ) %
      argv[3] % argv[2];

  std::cout <<
    boost::format( "  -----  -----  -----    ------------" );

  std::cout << std::endl;

  //============================================================================
  // Loop over species.
  //============================================================================

  nnt::species_list_t species_list = nnt::make_species_list( p_nuc );

  BOOST_FOREACH( nnt::Species species, species_list )
  {
    
    int i_z = Libnucnet__Species__getZ( species.getNucnetSpecies() );

    int i_a = Libnucnet__Species__getA( species.getNucnetSpecies() );

    int i_z2 = i_z - delta * Libnucnet__Species__getZ( p_sep );

    int i_a2 = i_a - delta * Libnucnet__Species__getA( p_sep );

    if( ( i_z2 > 0 && i_a2 > 0 ) or ( i_z2 == 0 && i_a2 == 1 ) )
    {
      p_species2 =
        Libnucnet__Nuc__getSpeciesByZA(
          p_nuc,
          i_z2,
          i_a2,
          ""
        );
    }
    else
    {
      p_species2 = NULL;
    }

    if( p_species2 )
    {

      std::cout <<
        boost::format( "%5d  %5d  %5d      %g\n" ) %
        i_z %
        (i_a - i_z) %
        i_a %
        (
          (double) delta * Libnucnet__Species__getMassExcess( p_sep ) +
          Libnucnet__Species__getMassExcess( p_species2 ) -
          Libnucnet__Species__getMassExcess( species.getNucnetSpecies() )
        );

    }
    else
    {
      if( B_OUTPUT_ALL_SPECIES )
      {

        std::cout <<
          boost::format( "%5d  %5d  %5d      %s\n" ) %
          i_z %
          (i_a - i_z) %
          i_a %
          S_BLANK;

      }

    }

  }

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__Nuc__free( p_nuc );

  return EXIT_SUCCESS;

}
