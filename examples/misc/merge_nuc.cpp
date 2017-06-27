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
//! \brief Example code to update nuclei from a first xml file from a second.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <fstream>
#include <iostream>

#include <nnt/auxiliary.h>

#include <Libnucnet__Nuc.h>

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
      " ../../data_pub/my_nuc.xml ../../data_pub/new_nuc.xml updated_nuc.xml" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc != 4 )
  {
    fprintf(
      stderr,
      "\nUsage: %s nuc1_xml nuc2_xml output_xml\n\n",
      argv[0]
    );
    fprintf(
      stderr, "  nuc1_xml = first nuclide xml filename\n\n"
    );
    fprintf(
      stderr, "  nuc2_xml = second nuclide xml filename\n\n"
    );
    fprintf(
      stderr, "  output_xml = output xml file\n\n"
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

  Libnucnet__Nuc * p_nuc;

  //============================================================================
  // Check input.
  //============================================================================

  check_input( argc, argv );

  //============================================================================
  // Read and store input.
  //============================================================================

  p_nuc = Libnucnet__Nuc__new_from_xml( argv[1], NULL );

  Libnucnet__Nuc__updateFromXml(
    p_nuc,
    argv[2],
    NULL
  );

  //============================================================================
  // Resort the species.
  //============================================================================

  Libnucnet__Nuc__setSpeciesCompareFunction(
    p_nuc,
    (Libnucnet__Species__compare_function) nnt::species_sort_by_z_then_a
  );

  Libnucnet__Nuc__sortSpecies( p_nuc );

  //============================================================================
  // Write to output.
  //============================================================================

  Libnucnet__Nuc__writeToXmlFile(
    p_nuc,
    argv[3]
  );

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__Nuc__free( p_nuc );
  return EXIT_SUCCESS;

}
