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
//       Example to demonstrate how to use libnucnet routines to parse
//       in a Libnucnet__Nuc structure of nuclear species from an input
//       xml file, compute and print out nuclear energies,
//       and clear the structure and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/
 
#include <Libnucnet__Nuc.h>

#define MY_BUFFER   32

/*##############################################################################
// Prototypes.
//############################################################################*/

void print_nuclei( Libnucnet__Nuc * );

int print_species( Libnucnet__Species *, Libnucnet__Nuc * );

double compute_species_mass( Libnucnet__Species * );

char *
compute_species_neutron_separation_energy(
  Libnucnet__Nuc *,
  Libnucnet__Species *
);

char *
compute_species_proton_separation_energy(
  Libnucnet__Nuc *,
  Libnucnet__Species *
);

/*##############################################################################
// main().
//############################################################################*/

int
main( int argc, char **argv ) {

  Libnucnet__Nuc *p_my_nuclei;

  if ( argc != 2 && argc != 3 ) {
      fprintf(
        stderr, "\nUsage: %s filename xpath\n\n", argv[0]
      );
      fprintf(
        stderr, "  filename = input nuclear data xml file or url.\n\n"
      );
      fprintf(
        stderr, "  xpath = xpath expression (optional).\n\n"
      );
      return EXIT_FAILURE;
   }

  /*============================================================================
  // Create nuclear species collection.
  //==========================================================================*/

  if( argc == 2 ) {
    p_my_nuclei = Libnucnet__Nuc__new_from_xml( argv[1], NULL );
  } else {
    p_my_nuclei = Libnucnet__Nuc__new_from_xml( argv[1], argv[2] );
  }
    
  /*============================================================================
  // Check that neutron and protons present.
  //==========================================================================*/
  
  if(
    !Libnucnet__Nuc__getSpeciesByName( p_my_nuclei, "h1" ) ||
    !Libnucnet__Nuc__getSpeciesByName( p_my_nuclei, "n" )
  )
  {
    fprintf( stderr, "Neutrons and protons must be in nuclear collection!\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Print the nuclear data.
  //==========================================================================*/

  print_nuclei( p_my_nuclei );

  /*============================================================================
  // Free the nuclear species collection.
  //==========================================================================*/

  Libnucnet__Nuc__free( p_my_nuclei );
  
  /*============================================================================
  // Done!
  //==========================================================================*/

  return EXIT_SUCCESS;

}

/*##############################################################################
// print_nuclei()
//############################################################################*/

void print_nuclei( Libnucnet__Nuc * self ) {

  /*============================================================================
  // Print out header.
  //==========================================================================*/

  fprintf(
    stdout,
    "\n  Z    A    Name  Mass (MeV)   B / A (MeV)   S_n (MeV)  S_p (MeV)\n"
  );
  fprintf(
    stdout,
    " ___  ___   ____  __________   ___________   _________  _________\n\n"
  );

  /*============================================================================
  // Print out species data.
  //==========================================================================*/

  Libnucnet__Nuc__iterateSpecies(
    self,
    (Libnucnet__Species__iterateFunction) print_species,
    self
  );

}

/*##############################################################################
// print_species()
//############################################################################*/

int
print_species(
  Libnucnet__Species *p_species,
  Libnucnet__Nuc *p_nuc
)
{

  double d_mass;
  char *s_sn, *s_sp;

  d_mass = compute_species_mass( p_species );

  s_sn =
    compute_species_neutron_separation_energy( p_nuc, p_species );

  s_sp =
    compute_species_proton_separation_energy( p_nuc, p_species );

  fprintf(
    stdout,
    "%4u %4u  %5s  %10.4f %10.4f       %s    %s\n",
    Libnucnet__Species__getZ( p_species ),
    Libnucnet__Species__getA( p_species ),
    Libnucnet__Species__getName( p_species ),
    d_mass,
    Libnucnet__Nuc__computeSpeciesBindingEnergy(
      p_nuc,
      p_species
    ) /
      Libnucnet__Species__getA( p_species ),
    s_sn,
    s_sp
  );

  free( s_sn );
  free( s_sp );

  return 1;

}

/*##############################################################################
// compute_species_mass().
//############################################################################*/

double compute_species_mass( Libnucnet__Species *p_species )
{

  return
    GSL_CONST_CGSM_UNIFIED_ATOMIC_MASS * 
    gsl_pow_2( GSL_CONST_CGSM_SPEED_OF_LIGHT ) *
    ( GSL_CONST_NUM_MICRO / GSL_CONST_CGSM_ELECTRON_VOLT ) *
    Libnucnet__Species__getA( p_species ) +
    Libnucnet__Species__getMassExcess( p_species );

}

/*##############################################################################
// compute_species_neutron_separation_energy().
//############################################################################*/

char *
compute_species_neutron_separation_energy(
  Libnucnet__Nuc *p_nuc,
  Libnucnet__Species *p_species
)
{

  Libnucnet__Species *p_species_2 = NULL;
  char *s_separation_energy;

  p_species_2 =
    Libnucnet__Nuc__getSpeciesByZA(
      p_nuc,
      Libnucnet__Species__getZ( p_species ),
      Libnucnet__Species__getA( p_species ) - 1,
      NULL
    );

  if( !p_species_2 )
    p_species_2 =
      Libnucnet__Nuc__getSpeciesByZA(
        p_nuc,
        Libnucnet__Species__getZ( p_species ),
        Libnucnet__Species__getA( p_species ) - 1,
        "g"
      );

  /*============================================================================
  // Return as string to allow for case when separation energy cannot be
  // computed.
  //==========================================================================*/

  s_separation_energy = ( char * ) malloc( sizeof( char ) * MY_BUFFER );

  if( p_species_2 )
    sprintf(
      s_separation_energy,
      "%.4f",
      Libnucnet__Nuc__computeSpeciesBindingEnergy( p_nuc, p_species ) -
        Libnucnet__Nuc__computeSpeciesBindingEnergy( p_nuc, p_species_2 )
    );
  else
    sprintf(
      s_separation_energy,
      "-------"
    );

  return s_separation_energy;

}
/*##############################################################################
// compute_species_proton_separation_energy().
//############################################################################*/

char *
compute_species_proton_separation_energy(
  Libnucnet__Nuc *p_nuc,
  Libnucnet__Species *p_species
)
{

  Libnucnet__Species *p_species_2;
  char *s_separation_energy;

  p_species_2 =
    Libnucnet__Nuc__getSpeciesByZA(
      p_nuc,
      Libnucnet__Species__getZ( p_species ) - 1,
      Libnucnet__Species__getA( p_species ) - 1,
      NULL
    );

  if( !p_species_2 )
    p_species_2 =
      Libnucnet__Nuc__getSpeciesByZA(
        p_nuc,
        Libnucnet__Species__getZ( p_species ),
        Libnucnet__Species__getA( p_species ) - 1,
        "g"
      );

  /*============================================================================
  // Return as string to allow for case when separation energy cannot be
  // computed.
  //==========================================================================*/

  s_separation_energy = ( char * ) malloc( sizeof( char ) * MY_BUFFER );

  if( p_species_2 )
    sprintf(
      s_separation_energy,
      "%.4f",
      Libnucnet__Nuc__computeSpeciesBindingEnergy( p_nuc, p_species ) -
        Libnucnet__Nuc__computeSpeciesBindingEnergy( p_nuc, p_species_2 )
    );
  else
    sprintf(
      s_separation_energy,
      "-------"
    );

  return s_separation_energy;

}
