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
//       equilibrium with a user-defined number density function and
//       integrand, and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/
 
#include <Libnuceq.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_const_cgsm.h>

/*##############################################################################
// defines.
//############################################################################*/

#define D_MIN_X_PRINT    1.e-10    /* Minimum mass fraction to print out. */

/*##############################################################################
// prototypes.
//############################################################################*/

int print_mass_fractions( Libnuceq__Species *, Libnuceq * );

double
electron_number_density_function(
  Libstatmech__Fermion *,
  double,
  double,
  void *
);

double
electron_number_density_integrand(
  Libstatmech__Fermion *,
  double,
  double,
  double,
  const char *
);

/*##############################################################################
// main().
//############################################################################*/

int main( int argc, char **argv ) {

  Libnucnet__Nuc *p_my_nuclei;
  Libnuceq *p_function, *p_integrand;

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
        stderr, "\n  shift: shift in chemical potential for integrand\n"
      );
      fprintf(
        stderr, "\n  nuc_xpath: XPath expression for nuclides\n\n"
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
  // Create the equilibria.
  //==========================================================================*/

  p_function = Libnuceq__new( p_my_nuclei );
  p_integrand  = Libnuceq__new( p_my_nuclei );

  /*============================================================================
  // Set the electron number density function and integrand.
  //==========================================================================*/

  Libnuceq__updateUserElectronNumberDensity(
    p_function,
    (Libstatmech__Fermion__Function) electron_number_density_function,
    NULL,
    NULL,
    NULL
  );

  Libnuceq__updateUserElectronNumberDensity(
    p_integrand,
    NULL,
    (Libstatmech__Fermion__Integrand) electron_number_density_integrand,
    NULL,
    argv[4]
  );

  /*============================================================================
  // Solve the equilibria.
  //==========================================================================*/

  Libnuceq__computeEquilibrium( p_function, atof( argv[2] ), atof( argv[3] ) );

  Libnuceq__computeEquilibrium( p_integrand, atof( argv[2] ), atof( argv[3] ) );

  /*============================================================================
  // Print abundances.
  //==========================================================================*/

  fprintf(
    stdout,
    "\nSpecies\t Mass fraction (function)\t Mass fraction (integrand)\n"
  );

  fprintf(
    stdout,
    "-------\t -----------------------\t --------------------------------\n\n"
  );

  Libnuceq__iterateSpecies(
    p_function,
    (Libnuceq__Species__iterateFunction) print_mass_fractions,
    p_integrand
  );

  /*============================================================================
  // Print basic parameters.
  //==========================================================================*/

  fprintf(
    stdout,
    "\nGlobal results for a T9 = %f, rho = %e g/cc WNSE:\n",
    Libnuceq__getT9( p_function ),
    Libnuceq__getRho( p_function )
  );

  /*============================================================================
  // Print function diagnostics.
  //==========================================================================*/

  fprintf(
    stdout,
    "\nDiagnostics for function calculation.\n"
  );

  fprintf(
    stdout,
    "\nYe = %f\n",
    Libnuceq__computeZMoment( p_function, 1 )
  );

  fprintf(
    stdout,
    "1 - Xsum = %e\n",
    1. - Libnuceq__computeAMoment( p_function, 1 )
  );

  fprintf(
    stdout,
    "\nNeutron chemical potential / kT = %f\n",
    Libnuceq__getMunkT( p_function )
  );

  fprintf(
    stdout,
    "Proton chemical potential / kT = %f\n",
    Libnuceq__getMupkT( p_function )
  );

  fprintf(
    stdout,
    "Electron chemical potential / kT = %f\n\n",
    Libnuceq__getMuekT( p_function )
  );

  /*============================================================================
  // Print diagnostics for integrand calculation.
  //==========================================================================*/

  fprintf(
    stdout,
    "\nDiagnostics for integrand calculation.\n"
  );

  fprintf(
    stdout,
    "\nYe = %f\n",
    Libnuceq__computeZMoment( p_integrand, 1 )
  );

  fprintf(
    stdout,
    "1 - Xsum = %e\n",
    1. - Libnuceq__computeAMoment( p_integrand, 1 )
  );

  fprintf(
    stdout,
    "\nNeutron chemical potential / kT = %f\n",
    Libnuceq__getMunkT( p_integrand )
  );

  fprintf(
    stdout,
    "Proton chemical potential / kT = %f\n",
    Libnuceq__getMupkT( p_integrand )
  );

  fprintf(
    stdout,
    "Electron chemical potential / kT = %f\n\n",
    Libnuceq__getMuekT( p_integrand )
  );

  /*============================================================================
  // Clean up.
  //==========================================================================*/

  Libnuceq__free( p_function );
  Libnuceq__free( p_integrand );
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
  Libnuceq *p_integrand
)
{

  double d_x, d_x_integrand;

  d_x =
    Libnuceq__Species__getAbundance( p_eq_species ) *
    Libnucnet__Species__getA(
      Libnuceq__Species__getNucSpecies( p_eq_species )
    );

  d_x_integrand =
    Libnuceq__Species__getAbundance(
      Libnuceq__getSpeciesByName(
        p_integrand,
        Libnucnet__Species__getName(
          Libnuceq__Species__getNucSpecies( p_eq_species )
        )
      )
    ) *
    Libnucnet__Species__getA(
      Libnuceq__Species__getNucSpecies( p_eq_species )
    );

  if( d_x > D_MIN_X_PRINT || d_x_integrand > D_MIN_X_PRINT )
    fprintf(
      stdout,
      "%s\t %e\t\t\t %e\n",
      Libnucnet__Species__getName(
        Libnuceq__Species__getNucSpecies( p_eq_species )
      ),
      d_x,
      d_x_integrand
    );

  return 1;

}

/*##############################################################################
// electron_number_density_integrand().
//############################################################################*/

double
electron_number_density_integrand(
  Libstatmech__Fermion *p_electron,
  double d_x,
  double d_T,
  double d_muekT,
  const char *s_shift
)
{

  return
    (
      Libstatmech__Fermion__getMultiplicity( p_electron ) /
      (
        sqrt( 2. ) *
        gsl_pow_2( M_PI ) *
        gsl_pow_3(
          GSL_CONST_CGSM_PLANCKS_CONSTANT_HBAR *
          GSL_CONST_CGSM_SPEED_OF_LIGHT
        )
      ) *
      pow(
        Libstatmech__Fermion__getRestMass( p_electron ) *
          GSL_CONST_CGSM_ELECTRON_VOLT *
          GSL_CONST_NUM_MEGA *
          GSL_CONST_CGSM_BOLTZMANN *
          d_T,
        3. / 2.
      )
    ) *
    sqrt( d_x ) / ( exp( d_x - d_muekT - atof( s_shift ) ) + 1. );

}

/*##############################################################################
// electron_number_density_function().
//############################################################################*/

double
electron_number_density_function(
  Libstatmech__Fermion *p_electron,
  double d_T,
  double d_muekT,
  void *p_data
)
{

  double d_result, d_fd;

  if( p_data ) {
    fprintf(
      stderr,
      "Invalid input!\n"
    );
    exit( EXIT_FAILURE );
  }

  if( d_muekT > -200 )
    d_fd = gsl_sf_fermi_dirac_half( d_muekT );
  else
    d_fd = exp( d_muekT );

  d_result =
    (
      Libstatmech__Fermion__getMultiplicity( p_electron ) /
      (
        sqrt( 2. ) *
        gsl_pow_2( M_PI ) *
        gsl_pow_3(
          GSL_CONST_CGSM_PLANCKS_CONSTANT_HBAR *
          GSL_CONST_CGSM_SPEED_OF_LIGHT
        )
      ) *
      pow(
        Libstatmech__Fermion__getRestMass( p_electron ) *
          GSL_CONST_CGSM_ELECTRON_VOLT *
          GSL_CONST_NUM_MEGA *
          GSL_CONST_CGSM_BOLTZMANN *
          d_T,
        3. / 2.
      )
    ) *
    gsl_sf_gamma( 1.5 ) * d_fd;

  return d_result;

}
