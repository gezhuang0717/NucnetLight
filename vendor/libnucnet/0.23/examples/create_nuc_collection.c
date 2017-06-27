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
//       Libnucnet__Nuc structure of nuclear species, add species to the
//       structure, remove species from the structure, retrieve data about
//       each species in the structure, and clear the structure and free the
//       allocated memory.
//     </abstract>
//   </description>
//
// </file>
//////////////////////////////////////////////////////////////////////////////*/
 
#include <Libnucnet__Nuc.h>

void print_nuclei( Libnucnet__Nuc * );

int print_species( Libnucnet__Species *, void * );

int
main( int argc, char **argv ) {

  Libnucnet__Nuc *p_my_nuclei;
  Libnucnet__Species *p_species, *p_copy1, *p_copy2;
  size_t i, i_partf_count = 10;
  unsigned i_z, i_a;
  int i_state;
  double d_mass_excess, d_spin;
  double a_t9[10] = {0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.};
  double a_log10_partf[10] = {1.,2.,3.,4.,5.,6.,7.,8.,9.,10.};
  gsl_vector *p_t9, *p_log10_partf;

  /*============================================================================
  // Check input.
  //==========================================================================*/

   if ( argc != 1 ) {
      fprintf(
        stderr, "\nUsage: %s\n\n", argv[0]
      );
      return EXIT_FAILURE;
   }

  /*============================================================================
  // Allocate and assign gsl vectors.
  //==========================================================================*/

  p_t9 = gsl_vector_alloc( i_partf_count );
  p_log10_partf = gsl_vector_alloc( i_partf_count );

  for( i = 0; i < i_partf_count; i++ ) {
    gsl_vector_set( p_t9, i, a_t9[i] );
    gsl_vector_set( p_log10_partf, i, a_log10_partf[i] );
  }

  /*============================================================================
  // Create collection of nuclear species.
  //==========================================================================*/

  p_my_nuclei = Libnucnet__Nuc__new( );

  print_nuclei( p_my_nuclei );

  /*============================================================================
  // Add neutrons.
  //==========================================================================*/

  fprintf( stdout, "\nAdd neutrons:\n\n" );

  i_z = 0; 
  i_a = 1;
  i_state = 0;
  d_mass_excess = 8.071;
  d_spin = 0.5;

  p_species =
    Libnucnet__Species__new(
      i_z,
      i_a,
      argv[0],
      i_state,
      NULL,
      d_mass_excess,
      d_spin,
      p_t9,
      p_log10_partf
  );

  if( !Libnucnet__Nuc__addSpecies( p_my_nuclei, p_species ) )
  {
    fprintf( stderr, "Couldn't add species.\n" );
    return EXIT_FAILURE;
  }

  print_nuclei( p_my_nuclei );

  /*============================================================================
  // Add di-neutrons.
  //==========================================================================*/

  fprintf( stdout, "\nAdd di-neutrons:\n\n" );

  i_z = 0; 
  i_a = 2;
  i_state = 0;
  d_mass_excess = 0.3;
  d_spin = 0.;

  p_species =
    Libnucnet__Species__new(
      i_z,
      i_a,
      "Made up data",
      i_state,
      NULL,
      d_mass_excess,
      d_spin,
      p_t9,
      p_log10_partf
  );

  if( !Libnucnet__Nuc__addSpecies( p_my_nuclei, p_species ) )
  {
    fprintf( stderr, "Couldn't add species.\n" );
    return EXIT_FAILURE;
  }

  print_nuclei( p_my_nuclei );

  /*============================================================================
  // Add protons.
  //==========================================================================*/

  fprintf( stdout, "\nAdd protons:\n\n" );

  i_z = 1; 
  i_a = 1;
  i_state = 0;
  d_mass_excess = 7.289;
  d_spin = 0.5;

  p_species =
    Libnucnet__Species__new(
      i_z,
      i_a,
      argv[0],
      i_state,
      NULL,
      d_mass_excess,
      d_spin,
      p_t9,
      p_log10_partf
    );

  if( !Libnucnet__Nuc__addSpecies( p_my_nuclei, p_species ) )
  {
    fprintf( stderr, "Couldn't add species.\n" );
    return EXIT_FAILURE;
  }

  /*========================================================================
  // Get two copies of protons.  One for an update.  One for later.
  //==========================================================================*/

  p_copy1 = Libnucnet__Species__copy( p_species );
  p_copy2 = Libnucnet__Species__copy( p_species );

  /*========================================================================
  // Update collection for illustration.
  //==========================================================================*/

  if( !Libnucnet__Nuc__updateSpecies( p_my_nuclei, p_copy1 ) )
  {
    fprintf( stderr, "Couldn't update species.\n" );
    return EXIT_FAILURE;
  }

  /*========================================================================
  // Print out data about the species.
  //==========================================================================*/

  print_nuclei( p_my_nuclei );

  /*============================================================================
  // Remove neutrons
  //==========================================================================*/

  fprintf( stdout, "\nRemove neutrons:\n\n" );

  Libnucnet__Nuc__removeSpecies(
    p_my_nuclei,
    Libnucnet__Nuc__getSpeciesByName( p_my_nuclei, "n" )
  );

  print_nuclei( p_my_nuclei );

  /*============================================================================
  // Free the collection.
  //==========================================================================*/

  printf( "\nFree the collection:\n\n" );

  Libnucnet__Nuc__free( p_my_nuclei );
  
  /*============================================================================
  // Recreate the collection.
  //==========================================================================*/

  p_my_nuclei = Libnucnet__Nuc__new();

  /*============================================================================
  // Add protons and print out data.
  //==========================================================================*/

  fprintf( stdout, "\nAdd protons again:\n\n" );

  if( !Libnucnet__Nuc__addSpecies( p_my_nuclei, p_copy2 ) )
  {
    fprintf( stderr, "Couldn't add species.\n" );
    return EXIT_FAILURE;
  }

  print_nuclei( p_my_nuclei );

  /*============================================================================
  // Final clean up
  //==========================================================================*/

  gsl_vector_free( p_t9 );
  gsl_vector_free( p_log10_partf );

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
  // Print out header information.
  //==========================================================================*/

  fprintf(
    stdout,
    "\n  Index  Z    A    Name    Mass Excess (MeV)    Spin   Data Source\n"
  );
  fprintf(
    stdout,
    "  _____ ___  ___  ______  ___________________   ____   __________\n\n"
  );

  /*============================================================================
  // Iterate to print species.
  //==========================================================================*/

  Libnucnet__Nuc__iterateSpecies(
    self,
    (Libnucnet__Species__iterateFunction) print_species,
    NULL
  );

  /*============================================================================
  // Print out the number of species.
  //==========================================================================*/

  fprintf(
    stdout,
    "\nThe collection has a total of %lu species.\n\n",
    (unsigned long) Libnucnet__Nuc__getNumberOfSpecies( self )
  );

  return;

}

/*##############################################################################
// print_species()
//############################################################################*/

int
print_species(
  Libnucnet__Species *p_species, void *p_user_data
)
{

  if( p_user_data ) {
    fprintf( stderr, "No user data should be passed to this function.\n" );
    return 0;
  }

  fprintf(
    stdout,
    "%5lu %4u %4u  %5s  %13.4f  %13.2f   %s\n",
    (unsigned long) Libnucnet__Species__getIndex( p_species ),
    Libnucnet__Species__getZ( p_species ),
    Libnucnet__Species__getA( p_species ),
    Libnucnet__Species__getName( p_species ),
    Libnucnet__Species__getMassExcess( p_species ),
    Libnucnet__Species__getSpin( p_species ),
    Libnucnet__Species__getSource( p_species )
  );

  return 1;

}

