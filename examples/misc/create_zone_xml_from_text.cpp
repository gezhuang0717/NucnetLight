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
//! \brief Example code to create a zone xml file from a text abundances file.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <Libnucnet.h>
#include "nnt/auxiliary.h"

//##############################################################################
// check_input().
//##############################################################################

void
check_input( int argc, char **argv )
{

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {
    std::cout << std::endl;
    std::cout << argv[0] <<
      " ../../data_pub/my_net.xml zone.txt ../../data_pub/zone.xml" <<
      " \"mass fractions\"" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if(
    argc < 5 ||
    argc > 6
  )
  {
    fprintf(
      stderr,
      "\nUsage: %s nuc_file zone_ascii zone_out flag xpath_nuc\n\n",
      argv[0]
    );
    fprintf(
      stderr, "  nuc_file = nuclear data xml filename\n\n"
    );
    fprintf(
      stderr, "  zone_ascii = input single zone ascii filename\n\n"
    );
    fprintf(
      stderr, "  zone_out = output zone xml filename\n\n"
    );
    fprintf(
      stderr, "  flag = zone ascii is either \"mass fractions\" or \"abundances\"\n\n"
    );
    fprintf(
      stderr,
      "  xpath_nuc = nuclear xpath expression (optional)\n\n"
    );
    std::cout << "For an example usage, type " << std::endl << std::endl;
    std::cout << argv[0] << " --example" << std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

}

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  double d_abundance;
  int i_z, i_a;
  Libnucnet *p_my_nucnet;
  nnt::Zone zone;
  Libnucnet__Species *p_species;
  FILE *p_file;

  //============================================================================
  // Check input.
  //============================================================================

  check_input( argc, argv );

  //============================================================================
  // Read and store input.
  //============================================================================

  p_my_nucnet = Libnucnet__new();

  if( argc ==5 )
    Libnucnet__Nuc__updateFromXml(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
      argv[1],
      NULL
    );  
  else
    Libnucnet__Nuc__updateFromXml(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
      argv[1],
      argv[5]
    );  

  //============================================================================
  // Create zone and add to my_nucnet. 
  //============================================================================

  zone.setNucnetZone(
    Libnucnet__Zone__new( 
      Libnucnet__getNet( p_my_nucnet ), 
      NULL,
      NULL,
      NULL
    ) 
  );

  Libnucnet__addZone( 
    p_my_nucnet,
    zone.getNucnetZone()
  );

  //============================================================================
  // Open file abundance file.
  //============================================================================

  p_file = fopen( argv[2], "r" );

  if( !p_file )
  {
    fprintf( stderr, "Couldn't open file %s.\n", argv[2] );
    return EXIT_FAILURE;
  }

  //============================================================================
  // Read in abundance data.
  //============================================================================

  do
  {
    fscanf( p_file, "%d %d %lf\n", &i_z, &i_a, &d_abundance );

    p_species =
      Libnucnet__Nuc__getSpeciesByZA(
        Libnucnet__Net__getNuc(
          Libnucnet__getNet( p_my_nucnet )
        ),
        i_z,
        i_a,
        NULL 
      );

    if( p_species )
    {

      if( strcmp( argv[4], "mass fractions" ) == 0 )
      {

        Libnucnet__Zone__updateSpeciesAbundance(
          zone.getNucnetZone(),
          p_species,
          d_abundance / i_a
        );

      }
      else if( strcmp( argv[4], "abundances" ) == 0 )
      {

        Libnucnet__Zone__updateSpeciesAbundance(
          zone.getNucnetZone(),
          p_species,
          d_abundance
        );

      }
      else
      {

        std::cerr << "Invalid abundance type." << std::endl;
        return EXIT_FAILURE;

      }

    }

  } while( !feof( p_file ) );

  fclose( p_file );

  //============================================================================
  // Normalize abundances. 
  //============================================================================

  zone.normalizeAbundances(); 

  //============================================================================
  // Write to output. 
  //============================================================================

  Libnucnet__updateZoneXmlMassFractionFormat(
    p_my_nucnet,
    "%.15e"
  );
        
  Libnucnet__writeZoneDataToXmlFile(
    p_my_nucnet,
    argv[3]
  );

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__free( p_my_nucnet );
  return EXIT_SUCCESS;

}
