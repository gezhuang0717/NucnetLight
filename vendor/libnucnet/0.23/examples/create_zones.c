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
//       a new Libnucnet structure from an
//       input xml file, print out data about the input zones,
//       remove the zones, add a zone, and clear the
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

int print_zone( Libnucnet__Zone *, int * );

int
print_species_abundance(
  Libnucnet__Species *, Libnucnet__Zone *
);

int zone_compare( const Libnucnet__Zone *, const Libnucnet__Zone * );

/*##############################################################################
// main().
//############################################################################*/

int main( int argc, char * argv[] ) {

  Libnucnet *p_my_nucnet;
  Libnucnet__Zone *p_zone, *p_copy;
  Libnucnet__Species *p_species;
  int i_count;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if ( argc!= 2 ) {
      fprintf(
        stderr, "\nUsage: %s filename\n\n", argv[0]
      );
      fprintf(
        stderr, "  filename = input nuclear network xml filename\n\n"
      );

      return EXIT_FAILURE;
  }

  /*============================================================================
  // Validate input file.
  //==========================================================================*/

  if( strcmp( VALIDATE, "yes" ) == 0 )
  {
    if( !Libnucnet__is_valid_input_xml( argv[1] ) ) {
      fprintf( stderr, "XML file is not valid Libnucnet input!\n" );
      return EXIT_FAILURE;
    }
  }

  /*============================================================================
  // Read input file.
  //==========================================================================*/

  p_my_nucnet = Libnucnet__new_from_xml( argv[1], NULL, NULL, NULL );

  /*============================================================================
  // Print number of zones.
  //==========================================================================*/

  
  fprintf(
    stdout,
    "\nNumber of zones = %lu\n\n",
    (unsigned long) Libnucnet__getNumberOfZones( p_my_nucnet )
  );

  /*============================================================================
  // Print zone abundances.  Set the zone_compare function to sort zones
  // alphabetically.
  //==========================================================================*/

  i_count = 1;

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function) zone_compare
  );

  Libnucnet__iterateZones(
    p_my_nucnet,
    (Libnucnet__Zone__iterateFunction) print_zone,
    &i_count
  );

  /*============================================================================
  // Remove zones.
  //==========================================================================*/

  fprintf( stdout, "Remove zones:\n\n" );

  Libnucnet__freeAllZones( p_my_nucnet );

  fprintf(
    stdout, "There are now %lu zones.\n\n",
    (unsigned long) Libnucnet__getNumberOfZones( p_my_nucnet )
  );

  /*============================================================================
  // Add a zone.
  //==========================================================================*/

  fprintf( stdout, "Now add a zone:\n" );

  p_zone =
    Libnucnet__Zone__new(
      Libnucnet__getNet( p_my_nucnet),
      "x3",
      "y3",
      "z3"
    );

  Libnucnet__addZone( p_my_nucnet, p_zone );

  /*============================================================================
  // Add a species abundance.
  //==========================================================================*/

  p_species =
    Libnucnet__Nuc__getSpeciesByName( 
      Libnucnet__Net__getNuc(
        Libnucnet__getNet( p_my_nucnet )
      ),
      "c12"
    );

  Libnucnet__Zone__updateSpeciesAbundance(
    p_zone, p_species, 1.0 / Libnucnet__Species__getA( p_species )
  );

  /*============================================================================
  // Print number of zones.
  //==========================================================================*/
  
  fprintf(
    stdout,
    "\nNumber of zones = %lu\n\n",
    (unsigned long) Libnucnet__getNumberOfZones( p_my_nucnet )
  );

  /*============================================================================
  // Print zone abundances.
  //==========================================================================*/

  i_count = 1;

  Libnucnet__iterateZones(
    p_my_nucnet,
    (Libnucnet__Zone__iterateFunction) print_zone,
    &i_count
  );

  /*============================================================================
  // Relabel the zone, change abundances, and print out the abundances again.
  //==========================================================================*/

  fprintf(
    stdout,
    "Relabel the zone, update abundances, and print out the abundances again:\n\n"
  );

  Libnucnet__relabelZone( p_my_nucnet, p_zone, "x4", "y4", NULL );

  p_species =
    Libnucnet__Nuc__getSpeciesByName( 
      Libnucnet__Net__getNuc(
        Libnucnet__getNet( p_my_nucnet )
      ),
      "c12"
    );

  Libnucnet__Zone__updateSpeciesAbundance(
    p_zone, p_species, 0.5 / Libnucnet__Species__getA( p_species )
  );

  p_species =
    Libnucnet__Nuc__getSpeciesByName( 
      Libnucnet__Net__getNuc(
        Libnucnet__getNet( p_my_nucnet )
      ),
      "o16"
    );

  Libnucnet__Zone__updateSpeciesAbundance(
    p_zone, p_species, 0.5 / Libnucnet__Species__getA( p_species )
  );

  i_count = 1;

  Libnucnet__iterateZones(
    p_my_nucnet,
    (Libnucnet__Zone__iterateFunction) print_zone,
    &i_count
  );

  /*============================================================================
  // Copy the zone, relabel the old zone, add the new zone, and print out.
  //==========================================================================*/

  fprintf(
    stdout,
    "Copy the zone, relabel the copied zone, and print out the abundances again:\n\n"
  );

  p_copy = Libnucnet__Zone__copy( p_zone );

  Libnucnet__relabelZone( p_my_nucnet, p_zone, "x5", "y5", NULL );

  Libnucnet__addZone( p_my_nucnet, p_copy );

  i_count = 1;

  Libnucnet__iterateZones(
    p_my_nucnet,
    (Libnucnet__Zone__iterateFunction) print_zone,
    &i_count
  );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}

/*##############################################################################
// print_zone()
//############################################################################*/

int
print_zone( Libnucnet__Zone *p_zone, int *p_count )
{

  fprintf(
    stdout,
    "     Zone %d",
    (*p_count)++
  );

  fprintf(
    stdout,
    " (Labels: %s %s %s)\n\n",
    Libnucnet__Zone__getLabel( p_zone, 1 ),
    Libnucnet__Zone__getLabel( p_zone, 2 ),
    Libnucnet__Zone__getLabel( p_zone, 3 )
  );

  fprintf( stdout, "Species\t Abundance\tMass Fraction\n" );
  fprintf( stdout, "=======\t============\t=============\n" );

  Libnucnet__Nuc__iterateSpecies(
    Libnucnet__Net__getNuc(
      Libnucnet__Zone__getNet( p_zone )
    ),
    (Libnucnet__Species__iterateFunction) print_species_abundance,
    p_zone
  );

  fprintf(
    stdout,
    "\nSum of mass fractions = %e\n",
    Libnucnet__Zone__computeAMoment( p_zone, 1 )
  );

  fprintf( stdout, "\n\n" );

  return 1;

}

/*##############################################################################
// print_species_abundance().
//############################################################################*/

int
print_species_abundance(
  Libnucnet__Species *p_species,
  Libnucnet__Zone *p_zone
)
{

  double d_abund;

  d_abund = 
     Libnucnet__Zone__getSpeciesAbundance( p_zone, p_species );

  if( d_abund > 0 ) {
     fprintf(
       stdout,
       "%s\t%e\t%e\n",
       Libnucnet__Species__getName( p_species ),
       d_abund,
       d_abund * Libnucnet__Species__getA( p_species )
     ); 
  }

  return 1;

}

/*##############################################################################
// zone_compare()
//############################################################################*/

int
zone_compare(
  const Libnucnet__Zone *p_zone1, const Libnucnet__Zone *p_zone2
)
{

  int i;
  
  if( 
     ( i =
         strcmp(
           Libnucnet__Zone__getLabel( p_zone1, 1 ),
           Libnucnet__Zone__getLabel( p_zone2, 1 )
         )
     )
  )
     return i;
          
  if( 
     ( i =
         strcmp(
           Libnucnet__Zone__getLabel( p_zone1, 2 ),
           Libnucnet__Zone__getLabel( p_zone2, 2 )
         )
     )
  )
     return i;
          
  return
    strcmp(
      Libnucnet__Zone__getLabel( p_zone1, 3 ),
      Libnucnet__Zone__getLabel( p_zone2, 3 )
    );

}
