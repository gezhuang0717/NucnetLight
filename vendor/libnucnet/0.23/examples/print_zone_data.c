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
//       a new Libnucnet structure from an input xml file, print out data
//       about the input zones, and clear the
//       structure and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libnucnet.h>

int print_zone_abundances( Libnucnet__Zone *,  int * );

int
zone_compare( 
  const Libnucnet__Zone *,
  const Libnucnet__Zone *
);

void my_func( const char *, const char *, const char *, const char *, void * );

int
print_species_abundance(
  Libnucnet__Species *, Libnucnet__Zone *p_zone
);

/*##############################################################################
// main().
//############################################################################*/

int main( int argc, char * argv[] ) {

  Libnucnet *p_my_nucnet;
  int i_count = 0;

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
  // Read input file.
  //==========================================================================*/

  p_my_nucnet = Libnucnet__new_from_xml( argv[1], NULL, NULL, NULL );

  /*============================================================================
  // Print number of zones.
  //==========================================================================*/

  printf(
    "\nNumber of zones = %lu\n\n",
    (unsigned long) Libnucnet__getNumberOfZones( p_my_nucnet )
  );

  /*============================================================================
  // Print zone abundances.  Iterate zones alphabetically by labels.
  //==========================================================================*/

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function) zone_compare
  );

  Libnucnet__iterateZones(
    p_my_nucnet,
    (Libnucnet__Zone__iterateFunction) print_zone_abundances,
    &i_count
  );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  Libnucnet__free( p_my_nucnet );
  return EXIT_SUCCESS;

}

/*##############################################################################
// print_zone_abundances().
//############################################################################*/

int
print_zone_abundances( Libnucnet__Zone *p_zone, int *p_count )
{

  printf(
    "Zone %d (Label 1 = %s  Label 2 = %s  Label 3 = %s)\n\n",
    (*p_count)++,
    Libnucnet__Zone__getLabel( p_zone, 1 ),
    Libnucnet__Zone__getLabel( p_zone, 2 ),
    Libnucnet__Zone__getLabel( p_zone, 3 )
  );
 
  printf( "Species\tAbundance\n" );
  printf( "=======\t=========\n" );

  Libnucnet__Nuc__iterateSpecies(
    Libnucnet__Net__getNuc(
      Libnucnet__Zone__getNet( p_zone )
    ),
    (Libnucnet__Species__iterateFunction) print_species_abundance,
     p_zone
  );

  printf( "\n\n" );

 return 1;

}

/*##############################################################################
// print_species_abundance().
//############################################################################*/

int
print_species_abundance
(
  Libnucnet__Species *p_species,
  Libnucnet__Zone *p_zone
)
{

  double d_abund;

  d_abund = 
     Libnucnet__Zone__getSpeciesAbundance( p_zone, p_species );

  if( d_abund > 0 )
    printf(
      "%s\t%f\n",
      Libnucnet__Species__getName( p_species ),
      d_abund
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
