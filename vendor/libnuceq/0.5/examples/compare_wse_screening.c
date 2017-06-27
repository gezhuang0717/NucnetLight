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
//       Example to demonstrate how to use libnuceq routines to create a
//       an equilibrium structure, compute weak nuclear statistical
//       equilibrium with and without a Coulomb correctio, and free the
//       allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/
 
#include <Libnuceq.h>
#include "coul_corr.h"

/*##############################################################################
// defines.
//############################################################################*/

#define D_MIN_X_PRINT    1.e-10    /* Minimum mass fraction to print out. */

/*##############################################################################
// prototypes.
//############################################################################*/

int print_mass_fractions( Libnuceq__Species *, Libnuceq * );

/*##############################################################################
// main().
//############################################################################*/

int main( int argc, char **argv ) {

  Libnucnet__Nuc *p_my_nuclei;
  Libnuceq *p_my_equil, *p_my_equil_corr;
  user_coul_corr_data my_corr_data = {-0.9052,0.6322};

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if ( argc < 4 || argc > 5) {
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
        stderr, "\n  nuc_xpath: XPath expression for nuclides (optional)\n\n"
      );
      return EXIT_FAILURE;
   }

  /*============================================================================
  // Get the nuclide collection.
  //==========================================================================*/

  if( argc == 5 )
    p_my_nuclei = Libnucnet__Nuc__new_from_xml( argv[1], argv[4] );
  else
    p_my_nuclei = Libnucnet__Nuc__new_from_xml( argv[1], NULL );

  /*============================================================================
  // Create equilibria.
  //==========================================================================*/

  p_my_equil = Libnuceq__new( p_my_nuclei );

  p_my_equil_corr = Libnuceq__new( p_my_nuclei );

  /*============================================================================
  // Set correction.
  //==========================================================================*/

  Libnuceq__setNseCorrectionFactorFunction(
    p_my_equil_corr,
    (Libnucnet__Species__nseCorrectionFactorFunction) my_coulomb_correction,
    &my_corr_data
  );

  /*============================================================================
  // Solve the equilibria.
  //==========================================================================*/

  Libnuceq__computeEquilibrium( p_my_equil, atof( argv[2] ), atof( argv[3] ) );

  Libnuceq__computeEquilibrium(
    p_my_equil_corr,
    atof( argv[2] ),
    atof( argv[3] )
  );

  /*============================================================================
  // Print abundances.
  //==========================================================================*/

  Libnuceq__iterateSpecies(
    p_my_equil,
    (Libnuceq__Species__iterateFunction) print_mass_fractions,
    p_my_equil_corr
  );

  /*============================================================================
  // Print out diagnostics.
  //==========================================================================*/

  fprintf(
    stdout,
    "\nGlobal results for a T9 = %f, rho = %e g/cc WNSE:\n",
    Libnuceq__getT9( p_my_equil ),
    Libnuceq__getRho( p_my_equil )
  );

  fprintf(
    stdout,
    "\nWithout screening:\n"
  );

  fprintf(
    stdout,
    "\nYe = %f\n",
    Libnuceq__computeZMoment( p_my_equil, 1 )
  );

  fprintf(
    stdout,
    "1 - Xsum = %e\n",
    1. - Libnuceq__computeAMoment( p_my_equil, 1 )
  );

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

  fprintf(
    stdout,
    "Electron chemical potential / kT = %f\n\n",
    Libnuceq__getMuekT( p_my_equil )
  );

  fprintf(
    stdout,
    "\nWith screening:\n"
  );

  fprintf(
    stdout,
    "\nYe = %f\n",
    Libnuceq__computeZMoment( p_my_equil_corr, 1 )
  );

  fprintf(
    stdout,
    "1 - Xsum = %e\n",
    1. - Libnuceq__computeAMoment( p_my_equil_corr, 1 )
  );

  fprintf(
    stdout,
    "\nNeutron chemical potential / kT = %f\n",
    Libnuceq__getMunkT( p_my_equil_corr )
  );

  fprintf(
    stdout,
    "Proton chemical potential / kT = %f\n",
    Libnuceq__getMupkT( p_my_equil_corr )
  );

  fprintf(
    stdout,
    "Electron chemical potential / kT = %f\n\n",
    Libnuceq__getMuekT( p_my_equil_corr )
  );

  /*============================================================================
  // Free equilibrium and nuclear collection.
  //==========================================================================*/

  Libnuceq__free( p_my_equil );
  Libnuceq__free( p_my_equil_corr );
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
print_mass_fractions(
  Libnuceq__Species *p_eq_species,
  Libnuceq *p_eq_corr
)
{

  double d_x, d_x_corr;
  Libnuceq__Species *p_eq_species_corr;

  p_eq_species_corr =
    Libnuceq__getSpeciesByName(
      p_eq_corr,
      Libnucnet__Species__getName(
        Libnuceq__Species__getNucSpecies( p_eq_species )
      )
    );

  d_x =
    Libnuceq__Species__getAbundance( p_eq_species ) *
    Libnucnet__Species__getA(
      Libnuceq__Species__getNucSpecies( p_eq_species )
    );

  d_x_corr =
    Libnuceq__Species__getAbundance( p_eq_species_corr ) *
    Libnucnet__Species__getA(
      Libnuceq__Species__getNucSpecies( p_eq_species_corr )
    );

  if( d_x > D_MIN_X_PRINT || d_x_corr > D_MIN_X_PRINT )
    fprintf(
      stdout,
      "Name: %s   Mass fraction: %10.4e  Mass fraction (Screened): %10.4e\n",
      Libnucnet__Species__getName(
        Libnuceq__Species__getNucSpecies( p_eq_species )
      ),
      d_x,
      d_x_corr
    );

  return 1;

}
