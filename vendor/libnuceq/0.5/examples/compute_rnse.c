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
//       Example to demonstrate how to use libnucnet routines to create a
//       restricted NSE.
//     </abstract>
//   </description>
// </file>
/////////////////////////////////////////////////////////////////////////////*/
 
#include <Libnuceq.h>
#include "coul_corr.h"

/*##############################################################################
// defines.
//############################################################################*/

#define D_MIN_X_PRINT  1.e-10
#define USE_CORRECTION "no"      /*  "yes" = use correction, "no" = don't */

/*##############################################################################
// prototypes.
//############################################################################*/

int print_mass_fractions( Libnuceq__Species *, Libnuceq * );

double rnse_prefactor( Libnuceq__Cluster *, Libnuceq__Species *, void *);

double rnse_constraint( Libnuceq__Species *, void *);

void mass_sum( Libnuceq__Species *, double * );

/*##############################################################################
// main().
//############################################################################*/

int main( int argc, char **argv ) {

  Libnucnet__Nuc *p_my_nuclei;
  Libnuceq *p_my_equil, *p_nse;
  Libnuceq__Cluster *p_cluster;
  double d_xsum = 0.;
  user_coul_corr_data my_corr_data = {-0.9052,0.6322};

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if ( argc < 6 || argc > 7 ) {
      fprintf(
        stderr, "\nUsage: %s input_file t9 rho Ye Yq nuc_xpath\n", argv[0]
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
        stderr, "\n  Xr: rNSE constraint\n"
      );
      fprintf(
        stderr, "\n  nuc_xpath: XPath expression for nuclides (optional)\n"
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

  if( atof( argv[4] ) < 0. || atof( argv[4] ) > 1. )
  {
    fprintf( stderr, "Invalid Ye!\n" );
    return EXIT_FAILURE;
  }

  if( atof( argv[5] ) < 0. || atof( argv[5] ) > 1. / 4. )
  {
    fprintf( stderr, "Invalid Xr!\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Create equilibrium.
  //==========================================================================*/

  p_my_equil = Libnuceq__new( p_my_nuclei );

  Libnuceq__setYe( p_my_equil, atof( argv[4] ) );

  p_cluster =
    Libnuceq__newCluster( p_my_equil, "[a > 3]" );

  Libnuceq__Cluster__updateConstraint( p_cluster, atof( argv[5] ) );

  Libnuceq__Cluster__updatePrefactorFunction(
    p_cluster,
    (Libnuceq__Cluster__prefactorFunction) rnse_prefactor,
    NULL
  );

  Libnuceq__Cluster__updateConstraintFunction(
    p_cluster,
    (Libnuceq__Cluster__constraint_function) rnse_constraint,
    NULL
  );

  /*============================================================================
  // For comparison, compute regular NSE.
  //==========================================================================*/
  
  p_nse = Libnuceq__new( p_my_nuclei );

  Libnuceq__setYe( p_nse, atof( argv[4] ) );

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
      p_nse,
      (Libnucnet__Species__nseCorrectionFactorFunction) my_coulomb_correction,
      &my_corr_data
    );
  }

  /*============================================================================
  // Compute the equilibria.
  //==========================================================================*/

  Libnuceq__computeEquilibrium( p_my_equil, atof( argv[2] ), atof( argv[3] ) );

  Libnuceq__computeEquilibrium( p_nse, atof( argv[2] ), atof( argv[3] ) );

  /*============================================================================
  // Print input parameters.
  //==========================================================================*/
  
  fprintf(
    stdout,
    "\nInput parameters: T9 = %f, rho = %e (g/cc), Ye = %f.\n",
    Libnuceq__getT9( p_my_equil ), 
    Libnuceq__getRho( p_my_equil ), 
    Libnuceq__computeZMoment( p_my_equil, 1 )
  );
 
  /*============================================================================
  // Print abundances.
  //==========================================================================*/
 
  fprintf(
    stdout,
    "\nSpecies    Mass Fraction     Nse Mass Fraction\n"
  );  

  fprintf(
    stdout,
    "-------    -------------     -----------------\n"
  );  

  Libnuceq__iterateSpecies(
    p_my_equil,
    (Libnuceq__Species__iterateFunction) print_mass_fractions,
    p_nse
  );

  /*============================================================================
  // Print diagnostics.
  //==========================================================================*/

  fprintf(
    stdout,
    "\nYe = %f\n",
    Libnuceq__computeZMoment( p_my_equil, 1 )
  );

  fprintf(
    stdout,
    "1 - xsum = %e\n",
    1. - Libnuceq__computeAMoment( p_my_equil, 1 )
  );

  Libnuceq__Cluster__iterateSpecies(
    p_cluster,
    (Libnuceq__Species__iterateFunction) mass_sum,
    &d_xsum
  );

  fprintf(
    stdout,
    "X (cluster) = %e\n",
    d_xsum
  );

  fprintf(
    stdout,
    "Rn = %e   Rp = %e\n",
    Libnuceq__Species__getAbundance(
      Libnuceq__getSpeciesByName( p_my_equil, "n" )
    ) /
    Libnuceq__Species__getAbundance(
      Libnuceq__getSpeciesByName( p_nse, "n" )
    ),
    Libnuceq__Species__getAbundance(
      Libnuceq__getSpeciesByName( p_my_equil, "h1" )
    ) /
    Libnuceq__Species__getAbundance(
      Libnuceq__getSpeciesByName( p_nse, "h1" )
    )
  );

  /*============================================================================
  // Print chemical potentials.
  //==========================================================================*/

  printf(
    "\nmu_n/kT = %e\nmu_p/kT = %e\nlambda = %e\n\n",
    Libnuceq__getMunkT( p_my_equil ),
    Libnuceq__getMupkT( p_my_equil ),
    Libnuceq__Cluster__getMukT( p_cluster )
  );

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  Libnuceq__free( p_my_equil );
  Libnuceq__free( p_nse );

  Libnucnet__Nuc__free( p_my_nuclei );

  /*============================================================================
  // Done!
  //==========================================================================*/

  return EXIT_SUCCESS;

}

/*##############################################################################
// rnse_prefactor()
//############################################################################*/

double
rnse_prefactor(
  Libnuceq__Cluster *p_cluster, Libnuceq__Species *p_species, void *p_data
)
{

  if( p_data )
  {
    fprintf(
      stderr,
      "No extra data to rnse_prefactor!\n"
    );
    exit( EXIT_FAILURE );
  }

  return
    -Libnuceq__Cluster__getMukT( p_cluster ) *
     Libnucnet__Species__getA(
       Libnuceq__Species__getNucSpecies( p_species )
     );

}

/*##############################################################################
// rnse_constraint()
//############################################################################*/

double
rnse_constraint(
  Libnuceq__Species *p_species, void *p_data
)
{

  if( p_data )
  {
    fprintf(
      stderr,
      "No extra data to rnse_prefactor!\n"
    );
    exit( EXIT_FAILURE );
  }

  return
    Libnuceq__Species__getAbundance( p_species ) *
    Libnucnet__Species__getA(
      Libnuceq__Species__getNucSpecies( p_species )
    );

}

/*##############################################################################
// mass_sum()
//############################################################################*/

void mass_sum( Libnuceq__Species *p_species, double *p_xsum )
{

  *p_xsum += rnse_constraint( p_species, NULL );

}

/*##############################################################################
// print_mass_fractions()
//############################################################################*/

int
print_mass_fractions( Libnuceq__Species *p_eq_species, Libnuceq *p_nse )
{

  Libnuceq__Species *p_nse_species;
  double d_x, d_x_nse;

  p_nse_species =
    Libnuceq__getSpeciesByName(
      p_nse,
      Libnucnet__Species__getName(
        Libnuceq__Species__getNucSpecies( p_eq_species )
      )
    );

  d_x =
    Libnuceq__Species__getAbundance( p_eq_species ) *
    Libnucnet__Species__getA(
      Libnuceq__Species__getNucSpecies( p_eq_species )
    );

  d_x_nse =
    Libnuceq__Species__getAbundance( p_nse_species ) *
    Libnucnet__Species__getA(
      Libnuceq__Species__getNucSpecies( p_nse_species )
    );

  if( d_x > D_MIN_X_PRINT || d_x_nse > D_MIN_X_PRINT )
    fprintf(
      stdout,
      "%5s   %15.6e   %15.6e\n",
      Libnucnet__Species__getName(
        Libnuceq__Species__getNucSpecies( p_eq_species )
      ),
      d_x,
      d_x_nse
    );

  return 1;

}
