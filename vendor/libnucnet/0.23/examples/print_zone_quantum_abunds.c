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
//       the quantum abundances and chemical potentials (in Maxwell-Boltzmann
//       statistics for the given zone, and clear the structure and free the
//       allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libnucnet.h>

#define T9   "t9"
#define RHO  "rho"

void print_out_zone( Libnucnet__Zone * );
int
print_species_quantum_abundance( Libnucnet__Species *, Libnucnet__Zone * );

int main( int argc, char * argv[] ) {

  Libnucnet *p_my_nucnet;
  Libnucnet__Zone *p_zone;

  if ( argc != 8 ) {
      fprintf(
        stderr,
        "\nUsage: %s nuc_file t9 rho mass_file label1 label2 label3\n\n",
        argv[0]
      );
      fprintf(
        stderr, "  nuc_file = input nuclear network data xml filename\n\n"
      );
      fprintf(
        stderr, "  mass_file = input mass fraction xml filename\n\n"
      );
      fprintf(
        stderr, "  t9 = input temperature in 10^9 K\n\n"
      );
      fprintf(
        stderr, "  rho = input mass density in g/cc\n\n"
      );
      fprintf(
        stderr, "  label1 = first zone label\n\n"
      );
      fprintf(
        stderr, "  label2 = second zone label\n\n"
      );
      fprintf(
        stderr, "  label3 = third zone label\n\n"
      );

      return EXIT_FAILURE;
  }

  /*============================================================================
  // Read input file.
  //==========================================================================*/

  p_my_nucnet = Libnucnet__new();

  Libnucnet__Nuc__updateFromXml(
    Libnucnet__Net__getNuc(
      Libnucnet__getNet( p_my_nucnet )
    ),
    argv[1],
    NULL
  );

  Libnucnet__assignZoneDataFromXml( p_my_nucnet, argv[2], NULL );

  /*============================================================================
  // Retrieve zone.
  //==========================================================================*/

  printf( "\n" );

  p_zone = Libnucnet__getZoneByLabels( 
    p_my_nucnet, argv[5], argv[6], argv[7] 
  );

  if( !p_zone ) {

    fprintf( stderr, "Zone not found!\n" );
    return EXIT_FAILURE;

  }

  /*============================================================================
  // Attach temperature and density to zone.
  //==========================================================================*/

  Libnucnet__Zone__updateProperty( p_zone, T9, NULL, NULL, argv[3] );

  Libnucnet__Zone__updateProperty( p_zone, RHO, NULL, NULL, argv[4] );

  /*============================================================================
  // Print out results.
  //==========================================================================*/

  printf(
    "T9 = %s  rho = %s\n",
    argv[3],
    argv[4]
  );

  printf( "For zone with labels %s %s %s:\n\n",
          Libnucnet__Zone__getLabel( p_zone, 1 ),
          Libnucnet__Zone__getLabel( p_zone, 2 ),
          Libnucnet__Zone__getLabel( p_zone, 3 )
  );

  print_out_zone( p_zone );

  printf( "\n" );

  Libnucnet__free( p_my_nucnet );
  return EXIT_SUCCESS;

}

/*##############################################################################
// print_out_zone()
//############################################################################*/

void print_out_zone( Libnucnet__Zone *p_zone ) {

  printf( "Species\t      Y     \t     Y_Q    \t   mu'/kT    \n" );
  printf( "=======\t============\t============\t=============\n" );

  Libnucnet__Nuc__iterateSpecies(
    Libnucnet__Net__getNuc(
      Libnucnet__Zone__getNet( p_zone )
    ),
    (Libnucnet__Species__iterateFunction) print_species_quantum_abundance,
    p_zone
  );

}

/*##############################################################################
// print_species_quantum_abundance().
//############################################################################*/

int
print_species_quantum_abundance(
  Libnucnet__Species *p_species,
  Libnucnet__Zone *p_zone
)
{

  double d_abund, d_yq;

  d_abund = 
     Libnucnet__Zone__getSpeciesAbundance( p_zone, p_species );

  if( d_abund > 0 )
  {
     d_yq =
       Libnucnet__Species__computeQuantumAbundance(
         p_species,
         atof( Libnucnet__Zone__getProperty( p_zone, T9, NULL, NULL ) ),
         atof( Libnucnet__Zone__getProperty( p_zone, RHO, NULL, NULL ) )
       );

     printf(
       "%s\t%e\t%e\t%e\n",
       Libnucnet__Species__getName( p_species ),
       d_abund,
       d_yq,
       log( d_abund / d_yq )
     ); 
  }

  return 1;

}

