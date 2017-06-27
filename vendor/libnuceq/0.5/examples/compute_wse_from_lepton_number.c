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
//       equilibrium from the fixed lepton number, and free the allocated
//       memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/
 
#include <Libnuceq.h>

/*##############################################################################
// defines.
//############################################################################*/

#define D_MIN_X_PRINT    1.e-10    /* Minimum mass fraction to print out. */

/*##############################################################################
// prototypes.
//############################################################################*/

int print_mass_fractions( Libnuceq__Species *, void * );

double neutrino_chemical_potential( Libnuceq *, double * );

double compute_neutrino_chemical_potential( Libnuceq *, double );

/*##############################################################################
// main().
//############################################################################*/

int main( int argc, char **argv ) {

  Libnucnet__Nuc *p_my_nuclei;
  Libnuceq *p_my_equil;
  double d_Y_L, d_mu_nue_kT, d_Ye;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if ( argc < 5 || argc > 6) {
      fprintf(
        stderr, "\nUsage: %s input_file t9 rho mu_nue_kT nuc_xpath\n", argv[0]
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
        stderr, "\n  Y_L: net lepton number per nucleon\n"
      );
      fprintf(
        stderr, "\n  nuc_xpath: XPath expression for nuclides (optional)\n\n"
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

  /*============================================================================
  // Set the electron neutrino chemical potential / kT.  Set the WSE
  // correction function.
  //==========================================================================*/

  d_Y_L = atof( argv[4] );

  Libnuceq__updateWseCorrectionFunction(
    p_my_equil,
    (Libnuceq__wseCorrectionFunction) neutrino_chemical_potential,
    &d_Y_L
  );

  /*============================================================================
  // Solve the equilibrium.
  //==========================================================================*/

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
  // Print equilibrium YL, Ye, and Ynue.
  //==========================================================================*/

  d_Ye = Libnuceq__computeZMoment( p_my_equil, 1 );

  fprintf(
    stdout,
    "\nGlobal results for a T9 = %f, rho = %e g/cc WNSE:\n",
    Libnuceq__getT9( p_my_equil ),
    Libnuceq__getRho( p_my_equil )
  );

  fprintf(
    stdout,
    "\nYL = %g\n",
    d_Y_L
  );

  fprintf(
    stdout,
    "Ye = %g\n",
    d_Ye
  );

  fprintf(
    stdout,
    "Ynue = %g\n",
    d_Y_L - d_Ye
  );

  /*============================================================================
  // Print 1 - sum of mass fractions for check of calculation.
  //==========================================================================*/

  fprintf(
    stdout,
    "1 - Xsum = %e\n",
    1. - Libnuceq__computeAMoment( p_my_equil, 1 )
  );

  /*============================================================================
  // Print chemical potentials.
  //==========================================================================*/

  fprintf(
    stdout,
    "\nChemical potentials (not including rest mass):\n"
  );

  d_mu_nue_kT =
    compute_neutrino_chemical_potential( p_my_equil, d_Y_L );

  fprintf(
    stdout,
    "\n  Neutron chemical potential / kT = %g\n",
    Libnuceq__getMunkT( p_my_equil )
  );

  fprintf(
    stdout,
    "  Proton chemical potential / kT = %g\n",
    Libnuceq__getMupkT( p_my_equil )
  );

  fprintf(
    stdout,
    "  Electron chemical potential / kT = %g\n",
    Libnuceq__getMuekT( p_my_equil )
  );

  fprintf(
    stdout,
    "  Neutrino chemical potential / kT = %g\n\n",
    d_mu_nue_kT
  );

  /*============================================================================
  // Check.
  //==========================================================================*/

  fprintf(
    stdout,
    "Check on WSE relation:\n"
  );

  fprintf(
    stdout,
    "\n  mu_p/kT + mu_e/kT - mu_n/kT - mu_nue/kT = %g\n\n",
    Libnuceq__getMupkT( p_my_equil ) +
    Libnuceq__getMuekT( p_my_equil ) -
    Libnuceq__getMunkT( p_my_equil ) -
    d_mu_nue_kT +
    (
      Libnucnet__Species__getMassExcess(
        Libnuceq__Species__getNucSpecies(
          Libnuceq__getSpeciesByName(
            p_my_equil,
            "h1"
          )
        )
      ) -
      Libnucnet__Species__getMassExcess(
        Libnuceq__Species__getNucSpecies(
          Libnuceq__getSpeciesByName(
            p_my_equil,
            "n"
          )
        )
      )
    ) /
    (
      ( GSL_CONST_NUM_KILO / GSL_CONST_CGSM_ELECTRON_VOLT ) *
      GSL_CONST_CGSM_BOLTZMANN *
      Libnuceq__getT9( p_my_equil )
    )
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
// neutrino_chemical_potential().
//############################################################################*/

double
neutrino_chemical_potential(
  Libnuceq *p_equil,
  double *p_y_lepton
)
{

  return
    compute_neutrino_chemical_potential(
      p_equil,
      *p_y_lepton
    );

}

/*##############################################################################
// compute_neutrino_chemical_potential().
//############################################################################*/

double
compute_neutrino_chemical_potential(
  Libnuceq *p_equil,
  double d_y_L
)
{

  Libstatmech__Fermion *p_neutrino;
  double d_mu_nue_kT, d_ye, d_number_density_neutrinos;

  d_ye = Libnuceq__computeZMoment( p_equil, 1 );

  /*============================================================================
  // Keep Ye bounded.
  //==========================================================================*/

  if( d_ye > 1. )
    d_ye = 1;
  else if( d_ye < -1. )
    d_ye = -1;

  d_number_density_neutrinos =
      Libnuceq__getRho( p_equil ) *
        GSL_CONST_NUM_AVOGADRO * ( d_y_L - d_ye );

  /*============================================================================
  // Create neutrino with zero rest mass.
  //==========================================================================*/

  p_neutrino =
    Libstatmech__Fermion__new(
      "nu_e",
      0.,
      1,
      0.
    );

  /*============================================================================
  // Compute chemical potential of neutrinos.
  //==========================================================================*/

  d_mu_nue_kT =
    Libstatmech__Fermion__computeChemicalPotential(
      p_neutrino,
      Libnuceq__getT9( p_equil ) * GSL_CONST_NUM_GIGA,
      d_number_density_neutrinos,
      NULL,
      NULL
    );

  /*============================================================================
  // Clean up and return.
  //==========================================================================*/

  Libstatmech__Fermion__free( p_neutrino );

  return d_mu_nue_kT;

}
