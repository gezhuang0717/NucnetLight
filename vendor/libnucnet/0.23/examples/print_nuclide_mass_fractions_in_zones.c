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
//       a new Libnucnet structure from an input xml file, print out 
//       abundances for the input species in all zones, and clear the
//       structure and free the allocated memory.
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

int print_zone_abundance( Libnucnet__Zone *, Libnucnet__Species * );

int
zone_compare( const Libnucnet__Zone *, const Libnucnet__Zone * );

/*##############################################################################
// main().
//############################################################################*/

int main( int argc, char * argv[] ) {

  Libnucnet *p_my_nucnet;
  Libnucnet__Species *p_species;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if ( argc!= 4 ) {
      fprintf(
        stderr, "\nUsage: %s nuc_file mass_file species\n\n", argv[0]
      );
      fprintf(
        stderr, "  nuc_file = input nuclear data xml filename\n\n"
      );
      fprintf(
        stderr, "  mass_file = input mass fraction xml filename\n\n"
      );
      fprintf(
        stderr, "  species = name of species\n\n"
      );

      return EXIT_FAILURE;
  }

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if( strcmp( VALIDATE, "yes" ) == 0 )
  {
    if( !Libnucnet__Nuc__is_valid_input_xml( argv[1] ) )
    {
      fprintf( stderr, "Invalid input nuclear data!\n" );
      return EXIT_FAILURE;
    }

    if( !Libnucnet__is_valid_zone_data_xml( argv[2] ) )
    {
      fprintf( stderr, "Invalid input zone data!\n" );
      return EXIT_FAILURE;
    }
  }

  /*============================================================================
  // Read input data.  For illustration, do this the hard way by reading
  // the nuclear data and input mass fractions separately.
  // The simpler thing to do would be to call:
  // 
  // p_my_nucnet = Libnucnet__new_from_xml( file, NULL, NULL );
  //
  // where file is an appropriate libnucnet input xml file.
  //==========================================================================*/

  p_my_nucnet = Libnucnet__new( );

  Libnucnet__Nuc__updateFromXml(
    Libnucnet__Net__getNuc(
      Libnucnet__getNet( p_my_nucnet )
    ),
    argv[1],
    NULL
  );

  Libnucnet__assignZoneDataFromXml( p_my_nucnet, argv[2], NULL );

  /*============================================================================
  // Retrieve desired species.
  //==========================================================================*/

  p_species =
    Libnucnet__Nuc__getSpeciesByName(
      Libnucnet__Net__getNuc(
        Libnucnet__getNet( p_my_nucnet )
      ),
      argv[3]
    );

  if( !p_species ) {
    fprintf( stderr, "Species not found!\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Print header.
  //==========================================================================*/

  fprintf( stdout, "\n" );

  fprintf( stdout, "\tZone\n" );

  fprintf( stdout, "\n" );

  fprintf(
    stdout,
    "Label 1\tLabel 2\tLabel3  \t\tMass Fraction of %s\n", argv[3]
  );
  fprintf( stdout, "-------\t-------\t--------\t\t---------------------\n" );
  fprintf( stdout, "\n" );

  /*============================================================================
  // Print abundance in zones.  Iterate the zones alphabetically by
  // last name then by first name.
  //==========================================================================*/

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function) zone_compare
  );

  Libnucnet__iterateZones(
    p_my_nucnet,
    (Libnucnet__Zone__iterateFunction) print_zone_abundance,
    p_species
  ); 

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}

/*##############################################################################
// print_zone_abundance()
//############################################################################*/

int
print_zone_abundance(
  Libnucnet__Zone *p_zone,
  Libnucnet__Species *p_species
)
{

  double d_abund;

  d_abund = 
    Libnucnet__Zone__getSpeciesAbundance( p_zone, p_species );

  fprintf(
    stdout,
    "%s\t%s\t%s \t\t%f\n",
    Libnucnet__Zone__getLabel( p_zone, 1 ),
    Libnucnet__Zone__getLabel( p_zone, 2 ),
    Libnucnet__Zone__getLabel( p_zone, 3 ),
    d_abund * Libnucnet__Species__getA( p_species )
  ); 

  return 1;

}

/*##############################################################################
// zone_compare()
//############################################################################*/

int
zone_compare(
  const Libnucnet__Zone *p_zone1,
  const Libnucnet__Zone *p_zone2
)
{

  int i;
  
  if( 
     ( i =
         strcmp(
           Libnucnet__Zone__getLabel( p_zone1, 3 ),
           Libnucnet__Zone__getLabel( p_zone2, 3 )
         )
     )
  )
     return i;
          
  if( 
     ( i =
         strcmp(
           Libnucnet__Zone__getLabel( p_zone1, 1 ),
           Libnucnet__Zone__getLabel( p_zone2, 1 )
         )
     )
  )
     return i;
          
  return
    strcmp(
      Libnucnet__Zone__getLabel( p_zone1, 2 ),
      Libnucnet__Zone__getLabel( p_zone2, 2 )
    );

}
