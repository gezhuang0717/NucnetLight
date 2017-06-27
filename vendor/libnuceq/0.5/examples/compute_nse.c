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
//       Example to demonstrate how to use libnuceq routines to create an
//       equilibrium, compute nuclear statistical equilibrium, print out
//       the results, and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/
 
#include <Libnuceq.h>

/*##############################################################################
// defines.
//############################################################################*/

#define D_MIN_X_PRINT  1.e-10   /*  Minimum mass fraction to print out. */

/*##############################################################################
// prototypes.
//############################################################################*/

int print_mass_fractions( Libnuceq__Species *, void * );
int cluster_abundance( Libnuceq__Species *, double * );

/*##############################################################################
// main().
//############################################################################*/

int main( int argc, char **argv ) {

  Libnucnet__Nuc *p_my_nuclei;
  Libnuceq *p_my_equil;
  Libnuceq__Cluster *p_cluster;
  double d_yh = 0;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if ( argc < 5 || argc > 6) {
      fprintf(
        stderr, "\nUsage: %s input_file t9 rho nuc_xpath\n", argv[0]
      );
      fprintf(
        stderr, "\n  input_file: name of input xml file\n"
      );
      fprintf(
        stderr, "\n  t9: Temperature in 10^9 K\n"
      );
      fprintf(
        stderr, "\n  rho: mass density in g/cc\n"
      );
      fprintf(
        stderr, "\n  ye: electron-to-nucleon ratio\n"
      );
      fprintf(
        stderr, "\n  nuc_xpath: XPath expression for nuclides (optional)\n"
      );
      return EXIT_FAILURE;
   }

  /*============================================================================
  // Get the nuclide collection.
  //==========================================================================*/

  if( argc == 6 )
    p_my_nuclei = Libnucnet__Nuc__new_from_xml( argv[1], argv[5] );
  else
    p_my_nuclei = Libnucnet__Nuc__new_from_xml( argv[1], NULL );

  /*============================================================================
  // Create equilibrium.
  //==========================================================================*/

  p_my_equil = Libnuceq__new( p_my_nuclei );

  Libnuceq__setYe( p_my_equil, atof( argv[4] ) );

  Libnuceq__computeEquilibrium( p_my_equil, atof( argv[2] ), atof( argv[3] ) );

  /*============================================================================
  // Print abundances.
  //==========================================================================*/

  Libnuceq__iterateSpecies(
    p_my_equil,
    (Libnuceq__Species__iterateFunction) print_mass_fractions,
    NULL
  );

  /*============================================================================
  // Print mass fraction diagnostic.
  //==========================================================================*/

  fprintf(
    stdout,
    "\nGlobal results for a T9 = %f, rho = %e g/cc, Ye = %f NSE:\n",
    Libnuceq__getT9( p_my_equil ),
    Libnuceq__getRho( p_my_equil ),
    Libnuceq__computeZMoment( p_my_equil, 1 )
  );

  fprintf(
    stdout,
    "1 - Xsum = %e\n",
    1. - Libnuceq__computeAMoment( p_my_equil, 1 )
  );

  /*============================================================================
  // Print number of heavy nuclei (Z>=6).  To get the number of heavy nuclei,
  // create a cluster of nuclei with Z >= 6. Since we don't recompute the
  // equilibrium, the cluster starts off with the original abundances from the
  // NSE; thus, summing over the cluster abundances gives the number of nuclei
  // with Z >= 6 in the NSE.
  //==========================================================================*/
 
  p_cluster = Libnuceq__newCluster( p_my_equil, "[z >= 6]" );

  Libnuceq__Cluster__iterateSpecies(
    p_cluster,
    (Libnuceq__Species__iterateFunction) cluster_abundance,
    &d_yh
  );

  fprintf(
    stdout,
    "Number of heavy nuclei per nucleon = %e\n",
    d_yh
  );

  /*============================================================================
  // Print chemical potentials.
  //==========================================================================*/

  fprintf(
    stdout,
    "\nNeutron chemical potential / kT = %f\n",
    Libnuceq__getMunkT( p_my_equil )
  );

  fprintf(
    stdout,
    "Proton chemical potential / kT = %f\n",
    Libnuceq__getMupkT( p_my_equil )
  );

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
/*##############################################################################
// cluster_abundance()
//############################################################################*/

int
cluster_abundance( Libnuceq__Species *p_species, double *p_yh )
{

  *p_yh += Libnuceq__Species__getAbundance( p_species );

  return 1;

}
