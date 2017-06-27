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
//       Example to demonstrate how to use libnucnet routines to compute
//       an (n,gamma)-(gamma,n) equilibrium.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/
 
#include <Libnuceq.h>

/*##############################################################################
// defines.
//############################################################################*/

#define BUF_SIZE  32
#define D_MIN_X_PRINT  1.e-10

/*##############################################################################
// prototypes.
//############################################################################*/

int print_mass_fractions( Libnuceq__Species *, void * );

/*##############################################################################
// main()
//############################################################################*/

int main( int argc, char **argv ) {

  FILE *p_file;
  Libnucnet__Nuc *p_my_nuclei;
  Libnuceq *p_my_equil;
  Libnuceq__Cluster *p_cluster;
  unsigned int i, i_z, i_z_max = 0;
  double d_yz;
  char s_xpath[BUF_SIZE];

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if ( argc < 6 || argc > 7 ) {
      fprintf(
        stderr, "\nUsage: %s input_file Y_Z file t9 rho Ye nuc_xpath\n", argv[0]
      );
      fprintf(
        stderr, "\n  input_file: name of input nuc xml file\n"
      );
      fprintf(
        stderr,
        "\n  Y_Z file: name of input text file with elemental abundances\n"
      );
      fprintf(
        stderr, "\n  t9: Temperature in 10^9 K\n"
      );
      fprintf(
        stderr, "\n  rho: mass density in g/cc\n"
      );
      fprintf(
        stderr, "\n  Ye: electron-to-nucleon ratio\n"
      );
      fprintf(
        stderr, "\n  nuc_xpath: XPath expression for nuclides (optional)\n\n"
      );
      return EXIT_FAILURE;
   }

  /*============================================================================
  // Get the nuclide collection.
  //==========================================================================*/

  if( argc == 7 )
    p_my_nuclei = Libnucnet__Nuc__new_from_xml( argv[1], argv[6] );
  else
    p_my_nuclei = Libnucnet__Nuc__new_from_xml( argv[1], NULL );

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if( atof( argv[5] ) < 0. || atof( argv[5] ) > 1. )
  {
    fprintf( stderr, "Invalid Ye!\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Create equilibrium.
  //==========================================================================*/

  p_my_equil = Libnuceq__new( p_my_nuclei );

  Libnuceq__setYe( p_my_equil, atof( argv[5] ) );

  /*============================================================================
  // Create (n,g)-(g,n) equilibrium clusters.
  //==========================================================================*/

  p_file = fopen( argv[2], "r" );

  if( !p_file )
  {
    fprintf( stderr, "\nCouldn't open file.\n" );
    return EXIT_FAILURE;
  }

  while( !feof( p_file ) )
  {
    fscanf( p_file, "%u%lf\n", &i_z, &d_yz );
    sprintf( s_xpath, "[z = %u]", i_z );
    p_cluster = Libnuceq__newCluster( p_my_equil, s_xpath );
    Libnuceq__Cluster__updateConstraint( p_cluster, d_yz );
    if( i_z > i_z_max ) i_z_max = i_z;
  }

  fclose( p_file );

  /*============================================================================
  // Compute equilibrium.
  //==========================================================================*/

  Libnuceq__computeEquilibrium( p_my_equil, atof( argv[3] ), atof( argv[4] ) );

  /*============================================================================
  // Print abundances.
  //==========================================================================*/

  Libnuceq__iterateSpecies(
    p_my_equil,
    (Libnuceq__Species__iterateFunction) print_mass_fractions,
    NULL
  );

  /*============================================================================
  // Print chemical potentials.
  //==========================================================================*/

  printf(
    "\nmu_n/kT = %e\nmu_p/kT = %e\n",
    Libnuceq__getMunkT( p_my_equil ),
    Libnuceq__getMupkT( p_my_equil )
  );

  for( i = 2; i <= i_z_max; i++ )
  {
    sprintf( s_xpath, "[z = %u]", i );
    p_cluster = Libnuceq__getCluster( p_my_equil, s_xpath );
    printf(
      "Z = %u  mu_h/kT = %e\n\n",
      i,
      Libnuceq__Cluster__getMukT( p_cluster )
    );
  }

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  Libnuceq__free( p_my_equil );
  Libnucnet__Nuc__free( p_my_nuclei );

  /*============================================================================
  // Done!
  //==========================================================================*/

  return EXIT_SUCCESS;

}

/*##############################################################################
// print_mass_fractions()
//############################################################################*/

int
print_mass_fractions( Libnuceq__Species *p_eq_species, void *p_data )
{

  double d_x;

  if( p_data ) exit( 0 );

  d_x =
    Libnuceq__Species__getAbundance( p_eq_species ) *
    Libnucnet__Species__getA(
      Libnuceq__Species__getNucSpecies( p_eq_species )
    );

  if( d_x > D_MIN_X_PRINT )
    fprintf(
      stdout,
      "Name: %s   Mass fraction: %e\n",
      Libnucnet__Species__getName(
        Libnuceq__Species__getNucSpecies( p_eq_species )
      ),
      d_x
    );

  return 1;

}
