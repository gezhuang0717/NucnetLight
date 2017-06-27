/*//////////////////////////////////////////////////////////////////////////////
// <file type="public">
//   <license>
//     This file was originally written by Tianhong Yu and Bradley S. Meyer.
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
//       Example to demonstrate how to use libstatmech routines to create and
//       store photons, compute thermodynamic quantities, and free
//       allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/

#include <Libstatmech.h>

double
my_number_density_function(
  Libstatmech__Boson *, double, double, void *
); 

double
my_pressure_function(
  Libstatmech__Boson *, double, double, void *
); 

double
my_energy_density_function(
  Libstatmech__Boson *, double, double, void *
); 

double
my_entropy_density_function(
  Libstatmech__Boson *, double, double, void *
); 

int main( int argc, char **argv )
{

  Libstatmech__Boson * p_photon;

  if( argc != 2 ) {
    fprintf(
      stderr,
      "\nUsage: %s temperature \n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  temperature = temperature in K\n\n"
    );
    return EXIT_FAILURE;
  } 

  p_photon =
    Libstatmech__Boson__new(
      "photon",
      0.,
      2,
      0.
   );

  fprintf(
    stdout,
    "The number density is %e 1/cm^3\n",
    Libstatmech__Boson__computeQuantity(
      p_photon,
      S_NUMBER_DENSITY,
      atof( argv[1] ),
      0.,
      NULL, 
      NULL
    )
  );

  fprintf(
    stdout,
    "The energy density is %e ergs/cm^3\n",
    Libstatmech__Boson__computeQuantity(
      p_photon,
      S_ENERGY_DENSITY,
      atof( argv[1] ),
      0.,
      NULL, 
      NULL 
    )
  );

  fprintf(
    stdout,
    "The pressure is %e dynes/cm^2\n",
    Libstatmech__Boson__computeQuantity(
      p_photon,
      S_PRESSURE,
      atof( argv[1] ),
      0.,
      NULL,
      NULL
    )
  );

  fprintf(
    stdout,
    "The entropy density is %e dynes/(K*cm^3)\n",
    Libstatmech__Boson__computeQuantity(
      p_photon,
      S_ENTROPY_DENSITY,
      atof( argv[1] ),
      0.,
      NULL, 
      NULL
    )
  );

  /*============================================================================
  // Compute quantities from function.  Do this by updating quantities.
  //==========================================================================*/

  if(
    !Libstatmech__Boson__updateQuantity(
      p_photon,
      S_NUMBER_DENSITY,
      (Libstatmech__Boson__Function) my_number_density_function,
      NULL
    )
  )
  {
    fprintf( stderr, "Couldn't update quantity!\n" );
    exit( EXIT_FAILURE );
  }

  if(
    !Libstatmech__Boson__updateQuantity(
      p_photon,
      S_PRESSURE,
      (Libstatmech__Boson__Function) my_pressure_function,
      NULL
    )
  )
  {
    fprintf( stderr, "Couldn't update quantity!\n" );
    exit( EXIT_FAILURE );
  }

  if(
    !Libstatmech__Boson__updateQuantity(
      p_photon,
      S_ENERGY_DENSITY,
      (Libstatmech__Boson__Function) my_energy_density_function,
      NULL
    )
  )
  {
    fprintf( stderr, "Couldn't update quantity!\n" );
    exit( EXIT_FAILURE );
  }

  if(
    !Libstatmech__Boson__updateQuantity(
      p_photon,
      S_ENTROPY_DENSITY,
      (Libstatmech__Boson__Function) my_entropy_density_function,
      NULL
    )
  )
  {
    fprintf( stderr, "Couldn't update quantity!\n" );
    exit( EXIT_FAILURE );
  }

  /*============================================================================
  // Compute quantities.
  //==========================================================================*/

  fprintf(
    stdout,
    "With our function the number density is %e 1/cm^3\n",
    Libstatmech__Boson__computeQuantity(
      p_photon,
      S_NUMBER_DENSITY,
      atof( argv[1] ),
      0.,
      NULL,
      NULL
    )
  );

  fprintf(
    stdout,
    "With our function the energy density is %e ergs/cm^3\n",
    Libstatmech__Boson__computeQuantity(
      p_photon,
      S_ENERGY_DENSITY,
      atof( argv[1] ),
      0.,
      NULL,
      NULL 
    )
  );

  fprintf(
    stdout,
    "With our function the pressure is %e dynes/cm^2\n",
    Libstatmech__Boson__computeQuantity(
      p_photon,
      S_PRESSURE,
      atof( argv[1] ),
      0.,
      NULL,
      NULL
    )
  );

  fprintf(
    stdout,
    "With our function the entropy density is %e dynes/(K*cm^3)\n",
    Libstatmech__Boson__computeQuantity(
      p_photon,
      S_ENTROPY_DENSITY,
      atof( argv[1] ),
      0.,
      NULL, 
      NULL
    )
  );

  /*============================================================================
  // Clean up and exit.
  //==========================================================================*/

  Libstatmech__Boson__free( p_photon );
  return EXIT_SUCCESS;

}

/*##############################################################################
// my_number_density_function()
//############################################################################*/

double
my_number_density_function(
  Libstatmech__Boson *p_photon,
  double d_T,
  double d_alpha,
  void *p_user_data
)
{

  double d_f;

  if( d_alpha > 0. || d_alpha < 0. ) 
  {
    fprintf( stderr, "Chemical potential must be zero.\n" );
    exit( 1 );
  }

  if( p_user_data )
  {
    fprintf( stderr, "Routine should have no extra user data.\n" );
    exit( 1 );
  }

  d_f =
    Libstatmech__Boson__getMultiplicity( p_photon ) *
    gsl_pow_3( GSL_CONST_CGSM_BOLTZMANN * d_T ) *
    gsl_sf_gamma( 3. ) *
    gsl_sf_zeta( 3. ) /
    (
      2. *
      gsl_pow_2( M_PI ) *
      gsl_pow_3(
        GSL_CONST_CGSM_PLANCKS_CONSTANT_HBAR *
        GSL_CONST_CGSM_SPEED_OF_LIGHT
      )
    );

  return d_f;

}
/*##############################################################################
// my_pressure_function()
//############################################################################*/

double
my_pressure_function(
  Libstatmech__Boson *p_photon,
  double d_T,
  double d_alpha,
  void *p_user_data
)
{

  if( !p_photon || p_user_data || d_alpha < 0. || d_alpha > 0. )
  {
    fprintf( stderr, "Invalid input.\n" );
    exit( 1 );
  }

  return
    4 *
    ( GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT /
      GSL_CONST_CGSM_SPEED_OF_LIGHT
    ) *
    gsl_pow_4( d_T ) /
    3.;

}

/*##############################################################################
// my_energy_density_function()
//############################################################################*/

double
my_energy_density_function(
  Libstatmech__Boson *p_photon,
  double d_T,
  double d_alpha,
  void *p_void
)
{

  double d_f;
  double d_a;

  p_void = &d_alpha;

  if( d_alpha > 0. || d_alpha < 0. || !p_void )
  {
    fprintf( stderr, "Invalid input.\n" );
    exit( EXIT_FAILURE );
  }

  if( !p_photon || d_T <= 0. )
  {
    fprintf( stderr, "Wrong input.\n" );
    exit( EXIT_FAILURE );
  }

  d_a = 
    4 *
    GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT /
    GSL_CONST_CGSM_SPEED_OF_LIGHT;

  d_f =
    d_a * 
    gsl_pow_4( d_T );

  return d_f;

}

/*##############################################################################
// my_entropy_density_function()
//############################################################################*/

double
my_entropy_density_function(
  Libstatmech__Boson *p_photon,
  double d_T,
  double d_alpha,
  void *p_void
)
{

  double d_f;
  double d_a;

  if( !p_photon || d_T <= 0. || p_void )
  {
    fprintf( stderr, "Wrong input.\n" );
    exit( EXIT_FAILURE );
  }

  if( d_alpha > 0. || d_alpha < 0. )
  {
    fprintf( stderr, "Chemical potential should be = 0.\n" );
    exit( EXIT_FAILURE );
  }

  d_a = 
    4 *
    GSL_CONST_CGSM_STEFAN_BOLTZMANN_CONSTANT /
    GSL_CONST_CGSM_SPEED_OF_LIGHT;

  d_f =
    4. * 
    d_a * 
    gsl_pow_3( d_T ) /
    3.;

  return d_f;

}
