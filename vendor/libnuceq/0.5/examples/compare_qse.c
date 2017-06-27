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
//       Example to demonstrate how to use libnuceq routines to create two
//       equilibria, compute the quasi-statistical equilibrium for each,
//       print out the results, and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/
 
#include <Libnuceq.h>
#include "coul_corr.h"

/*##############################################################################
// defines.
//############################################################################*/

#define D_MIN_X_PRINT  1.e-10    /*  Minimum mass fraction to print out. */
#define USE_CORRECTION "no"      /*  "yes" = use correction, "no" = don't */

#define S_NUC_XPATH    "[z >= 6]"

/*##############################################################################
// prototypes.
//############################################################################*/

int print_mass_fractions( Libnuceq__Species *, void * );

void print_diagnostics( Libnuceq * );

/*##############################################################################
// main().
//############################################################################*/

int main( int argc, char **argv ) {

  Libnucnet__Nuc *p_my_nuclei;
  Libnuceq *p_my_equil, *p_my_equil2;
  Libnuceq__Cluster *p_cluster;
  user_coul_corr_data my_corr_data = {-0.9052,0.6322};

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if ( argc < 7 || argc > 8 ) {
      fprintf(
        stderr, "\nUsage: %s input_file t9 rho Ye Yq Yq2 nuc_xpath\n", argv[0]
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
        stderr, "\n  Ye: electron-to-nucleon ratio\n"
      );
      fprintf(
        stderr, "\n  Yq: QSE constraint\n"
      );
      fprintf(
        stderr, "\n  Yq2: QSE constraint 2\n"
      );
      fprintf(
        stderr, "\n  nuc_xpath: XPath expression for nuclides (optional)\n"
      );
      return EXIT_FAILURE;
   }

  /*============================================================================
  // Get the nuclide collection.
  //==========================================================================*/

  if( argc == 8 )
    p_my_nuclei = Libnucnet__Nuc__new_from_xml( argv[1], argv[7] );
  else
    p_my_nuclei = Libnucnet__Nuc__new_from_xml( argv[1], NULL );

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if( atof( argv[4] ) < 0. || atof( argv[4] ) > 1. )
  {
    fprintf( stderr, "Invalid Ye!\n" );
    return EXIT_FAILURE;
  }

  if(
    atof( argv[5] ) < 0. ||
    atof( argv[5] ) > 1. / 12.
  )
  {
    fprintf( stderr, "Invalid Yq!\n" );
    return EXIT_FAILURE;
  }

  if(
    atof( argv[6] ) < 0. ||
    atof( argv[6] ) > 1. / 12.
  )
  {
    fprintf( stderr, "Invalid Yq2!\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Create first equilibrium.
  //==========================================================================*/

  p_my_equil = Libnuceq__new( p_my_nuclei );

  Libnuceq__setYe( p_my_equil, atof( argv[4] ) );

  p_cluster = Libnuceq__newCluster( p_my_equil, S_NUC_XPATH );
  Libnuceq__Cluster__updateConstraint( p_cluster, atof( argv[5] ) ); 

  /*============================================================================
  // Create second equilibrium.  Simply copy cluster instead of creating new.
  //==========================================================================*/

  p_my_equil2 = Libnuceq__new( p_my_nuclei );

  Libnuceq__setYe( p_my_equil2, atof( argv[4] ) );

  Libnuceq__copy_clusters( p_my_equil2, p_my_equil );

  p_cluster = Libnuceq__getCluster( p_my_equil2, S_NUC_XPATH );
  Libnuceq__Cluster__updateConstraint( p_cluster, atof( argv[6] ) ); 

  /*============================================================================
  // Set correction if desired.
  //==========================================================================*/

  if( strcmp( USE_CORRECTION, "yes" ) )
  {
    Libnuceq__setNseCorrectionFactorFunction(
      p_my_equil,
      (Libnucnet__Species__nseCorrectionFactorFunction) my_coulomb_correction,
      &my_corr_data
    );
    Libnuceq__setNseCorrectionFactorFunction(
      p_my_equil2,
      (Libnucnet__Species__nseCorrectionFactorFunction) my_coulomb_correction,
      &my_corr_data
    );
  }

  /*============================================================================
  // Compute the equilibria.
  //==========================================================================*/

  Libnuceq__computeEquilibrium( p_my_equil, atof( argv[2] ), atof( argv[3] ) );

  Libnuceq__computeEquilibrium( p_my_equil2, atof( argv[2] ), atof( argv[3] ) );

  /*============================================================================
  // Print abundances.
  //==========================================================================*/

  fprintf( stdout, "  Species     X(Cluster 1)   X(Cluster 2)\n" );

  fprintf( stdout, "-----------   ------------   ------------\n" );

  Libnuceq__iterateSpecies(
    p_my_equil,
    (Libnuceq__Species__iterateFunction) print_mass_fractions,
    p_my_equil2
  );

  /*============================================================================
  // Print chemical potentials.
  //==========================================================================*/

  print_diagnostics( p_my_equil );

  print_diagnostics( p_my_equil2 );

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  Libnuceq__free( p_my_equil );
  Libnuceq__free( p_my_equil2 );
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

  double d_x, d_x2;

  Libnuceq * p_equil2 = (Libnuceq *) p_data;

  Libnuceq__Species * p_eq_species2 =
    Libnuceq__getSpeciesByName(
      p_equil2,
      Libnucnet__Species__getName(
        Libnuceq__Species__getNucSpecies( p_eq_species )
      )
    );

  d_x =
    Libnuceq__Species__getAbundance( p_eq_species ) *
    Libnucnet__Species__getA(
      Libnuceq__Species__getNucSpecies( p_eq_species )
    );

  d_x2 =
    Libnuceq__Species__getAbundance( p_eq_species2 ) *
    Libnucnet__Species__getA(
      Libnuceq__Species__getNucSpecies( p_eq_species2 )
    );

  if( d_x > D_MIN_X_PRINT || d_x2 > D_MIN_X_PRINT )
    fprintf(
      stdout,
      "Name: %5s   %e   %e\n",
      Libnucnet__Species__getName(
        Libnuceq__Species__getNucSpecies( p_eq_species )
      ),
      d_x,
      d_x2
    );

  return 1;

}

/*##############################################################################
// print_diagnostics()
//############################################################################*/

void
print_diagnostics( Libnuceq * p_equil )
{

  Libnuceq__Cluster * p_cluster =
    Libnuceq__getCluster( p_equil, S_NUC_XPATH );

  fprintf(
    stdout,
    "\nQSE for T9 = %f, rho = %e, Ye = %f,\n",
    Libnuceq__getT9( p_equil ),
    Libnuceq__getRho( p_equil ),
    Libnuceq__computeZMoment( p_equil, 1 )
  );

  fprintf(
    stdout,
    "and for a QSE cluster defined by XPath expression %s\n",
    Libnuceq__Cluster__getXPathString( p_cluster )
  );

  fprintf(
    stdout,
    "with abundance constraint Yh = %e:\n",
    Libnuceq__Cluster__getConstraint( p_cluster )
  );

  fprintf(
    stdout,
    "\nmu_n/kT = %e\nmu_p/kT = %e\nmu_h/kT = %e\n\n",
    Libnuceq__getMunkT( p_equil ),
    Libnuceq__getMupkT( p_equil ),
    Libnuceq__Cluster__getMukT( p_cluster )
  );

}
