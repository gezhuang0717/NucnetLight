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
//       a (p,gamma)-(gamma,p) equilibrium.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/
 
#include <Libnuceq.h>

#define BUF_SIZE  32

int print_abunds( Libnuceq__Species *, void * );

int
sort_by_n_function(
  const Libnucnet__Species *,
  const Libnucnet__Species *
);

int main( int argc, char **argv ) {

  FILE *p_file;
  Libnucnet__Nuc *p_my_nuclei;
  Libnuceq *p_my_equil;
  Libnuceq__Cluster *p_cluster;
  unsigned int i, i_n, i_n_max = 0;
  double d_yn;
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
        "\n  Y_N file: name of input text file with isotonic abundances\n"
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
        stderr, "\n  nuc_xpath: XPath expression for nuclides\n\n"
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
  // Create (p,g)-(g,p) equilibrium clusters.
  //==========================================================================*/

  p_file = fopen( argv[2], "r" );

  if( !p_file )
  {
    fprintf( stderr, "\nCouldn't open file.\n" );
    return EXIT_FAILURE;
  }

  while( !feof( p_file ) )
  {
    fscanf( p_file, "%u%lf\n", &i_n, &d_yn );
    sprintf( s_xpath, "[(a - z) = %u]", i_n );
    p_cluster = Libnuceq__newCluster( p_my_equil, s_xpath );
    Libnuceq__Cluster__updateConstraint( p_cluster, d_yn );
    if( i_n > i_n_max ) i_n_max = i_n;
  }

  fclose( p_file );

  /*============================================================================
  // Compute equilibrium.
  //==========================================================================*/

  Libnuceq__computeEquilibrium( p_my_equil, atof( argv[3] ), atof( argv[4] ) );

  /*============================================================================
  // Print abundances.
  //==========================================================================*/

  Libnucnet__Nuc__setSpeciesCompareFunction(
    Libnuceq__getNuc( p_my_equil ),
    (Libnucnet__Species__compare_function) sort_by_n_function
  );

  Libnucnet__Nuc__sortSpecies( Libnuceq__getNuc( p_my_equil ) );

  Libnuceq__iterateSpecies(
    p_my_equil,
    (Libnuceq__Species__iterateFunction) print_abunds,
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

  for( i = 2; i <= i_n_max; i++ )
  {
    sprintf( s_xpath, "[(a - z) = %u]", i );
    p_cluster = Libnuceq__getCluster( p_my_equil, s_xpath );
    printf(
      "N = %u  mu_h/kT = %e\n\n",
      i,
      Libnuceq__Cluster__getMukT( p_cluster )
    );
  }

  /*============================================================================
  // Free equilibrium.
  //==========================================================================*/

  Libnuceq__free( p_my_equil );

  /*============================================================================
  // Done!
  //==========================================================================*/

  return EXIT_SUCCESS;

}

/*##############################################################################
// print_abunds()
//############################################################################*/

int
print_abunds( Libnuceq__Species *p_eq_species, void *p_data )
{

  if( p_data ) exit( 0 );

  fprintf(
    stdout,
    "%s\t%u\t%u\t%u\t%e\n",
    Libnucnet__Species__getName(
      Libnuceq__Species__getNucSpecies( p_eq_species )
    ),
    Libnucnet__Species__getZ(
      Libnuceq__Species__getNucSpecies( p_eq_species )
    ),
    Libnucnet__Species__getA(
      Libnuceq__Species__getNucSpecies( p_eq_species )
    ) -
    Libnucnet__Species__getZ(
      Libnuceq__Species__getNucSpecies( p_eq_species )
    ),
    Libnucnet__Species__getA(
      Libnuceq__Species__getNucSpecies( p_eq_species )
    ),
    Libnuceq__Species__getAbundance( p_eq_species ) *
    Libnucnet__Species__getA(
      Libnuceq__Species__getNucSpecies( p_eq_species )
    )
  );

  return 1;

}

/*##############################################################################
// sort_by_n_function()
//############################################################################*/

int
sort_by_n_function(
  const Libnucnet__Species *p_species1,
  const Libnucnet__Species *p_species2
)
{

  int i;

  if( 
      Libnucnet__Species__getA( p_species1 ) -
        Libnucnet__Species__getZ( p_species1 ) < 
      Libnucnet__Species__getA( p_species2 ) -
        Libnucnet__Species__getZ( p_species2 )
  )
    return -1;
  else if(
      Libnucnet__Species__getA( p_species1 ) -
        Libnucnet__Species__getZ( p_species1 ) >
      Libnucnet__Species__getA( p_species2 ) -
        Libnucnet__Species__getZ( p_species2 )
  )
    return 1;

  if( 
      Libnucnet__Species__getA( p_species1 ) <
      Libnucnet__Species__getA( p_species2 )
  )
    return -1;
  else if(
      Libnucnet__Species__getA( p_species1 ) > 
      Libnucnet__Species__getA( p_species2 )
  )
    return 1;

  i =
    strcmp(
      Libnucnet__Species__getName( p_species1 ),
      Libnucnet__Species__getName( p_species2 )
    );

  if( i == 0 ) {
    return 0;
  } else {
    return -GSL_SIGN( i );
  }

}
   
