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
//       structure, write the data to an xml file, and clear the structure
//       and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/
 
#include <Libnuceq.h>

int print_abunds( Libnuceq__Species *, void * );

int main( int argc, char **argv ) {

  Libnucnet__Nuc *p_my_nuclei;
  Libnuceq *p_my_equil;
  Libnuceq__Cluster *p_cluster;

  /*============================================================================
  // Check input.
  //==========================================================================*/

  if ( argc < 5 || argc > 6 ) {
      fprintf(
        stderr, "\nUsage: %s input_file t9 rho Yq nuc_xpath\n", argv[0]
      );
      fprintf(
        stderr, "\n  input_file: name of input text file\n"
      );
      fprintf(
        stderr, "\n  t9: Temperature in 10^9 K\n"
      );
      fprintf(
        stderr, "\n  rho: mass density in g/cc\n"
      );
      fprintf(
        stderr, "\n  Yq: QSE constraint\n"
      );
      fprintf(
        stderr, "\n  nuc_xpath: XPath expression for nuclides\n"
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
  // Check input.
  //==========================================================================*/

  if( atof( argv[4] ) < 0. || atof( argv[4] ) > 1. / 12. )
  {
    fprintf( stderr, "Invalid Yq!\n" );
    return EXIT_FAILURE;
  }

  /*============================================================================
  // Create equilibrium.
  //==========================================================================*/

  p_my_equil = Libnuceq__new( p_my_nuclei );

  p_cluster = Libnuceq__newCluster( p_my_equil, "[z >= 6]" );
  Libnuceq__Cluster__updateConstraint( p_cluster, atof( argv[4] ) );

  Libnuceq__computeEquilibrium( p_my_equil, atof( argv[2] ), atof( argv[3] ) );

  /*============================================================================
  // Print abundances.
  //==========================================================================*/

  Libnuceq__iterateSpecies(
    p_my_equil,
    (Libnuceq__Species__iterateFunction) print_abunds,
    NULL
  );

  /*============================================================================
  // Print equilibrium Ye.
  //==========================================================================*/

  fprintf(
    stdout,
    "Ye = %f\n",
    Libnuceq__computeZMoment( p_my_equil, 1 )
  );

  /*============================================================================
  // Print chemical potentials.
  //==========================================================================*/

  p_cluster = Libnuceq__getCluster( p_my_equil, "[z >= 6]" );

  fprintf(
    stdout,
    "\nmu_n/kT = %e\nmu_p/kT = %e\nmu_h/kT = %e\n\n",
    Libnuceq__getMunkT( p_my_equil ),
    Libnuceq__getMupkT( p_my_equil ),
    Libnuceq__Cluster__getMukT( p_cluster )
  );

  /*============================================================================
  // Free equilibrium.
  //==========================================================================*/

  Libnuceq__free( p_my_equil );

  /*============================================================================
  // Done!
  //==========================================================================*/

  return EXIT_SUCCESS;

}

int
print_abunds( Libnuceq__Species *p_eq_species, void *p_data )
{

  if( p_data ) exit( 0 );

  fprintf(
    stdout,
    "Name: %s   Mass fraction: %e\n",
    Libnucnet__Species__getName(
      Libnuceq__Species__getNucSpecies( p_eq_species )
    ),
    p_eq_species->dAbundance * p_eq_species->pNucSpecies->iA
  );

  return 1;

}
