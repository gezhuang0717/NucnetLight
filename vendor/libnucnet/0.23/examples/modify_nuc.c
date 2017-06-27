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
//       xml file, create a view and modify data included in the view,
//       and clear the structures and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/
 
#include <Libnucnet__Nuc.h>

int modify_species( Libnucnet__Species *, double * );

void print_nuclei( Libnucnet__Nuc * );

int print_species( Libnucnet__Species *, void * );

int
main( int argc, char **argv ) {

  Libnucnet__Nuc *p_my_nuclei;
  Libnucnet__NucView *p_view;
  double d_addend;

  if ( argc != 4 ) {
      fprintf(
        stderr, "\nUsage: %s filename xpath addend\n\n", argv[0]
      );
      fprintf(
        stderr, "  filename = input nuclear data xml file or url.\n\n"
      );
      fprintf(
        stderr, "  xpath = xpath expression for view nuclides.\n\n"
      );
      fprintf(
        stderr, "  addend = addend to mass excess.\n\n"
      );
      return EXIT_FAILURE;
   }

  /*============================================================================
  // Create nuclear species collection.
  //==========================================================================*/

  p_my_nuclei = Libnucnet__Nuc__new_from_xml( argv[1], NULL );
    
  /*============================================================================
  // Get view.
  //==========================================================================*/

  p_view = Libnucnet__NucView__new( p_my_nuclei, argv[2] );

  /*============================================================================
  // Modify source info in view nuclides.
  //==========================================================================*/

  d_addend = atof( argv[3] );

  Libnucnet__Nuc__iterateSpecies(
    Libnucnet__NucView__getNuc( p_view ),
    (Libnucnet__Species__iterateFunction) modify_species,
    &d_addend
  );

  /*============================================================================
  // Print the nuclear data.
  //==========================================================================*/

  print_nuclei( p_my_nuclei );

  /*============================================================================
  // Print the number of nuclei whose data were modified.
  //==========================================================================*/

  fprintf(
    stdout,
    "The number of nuclei with modified data (those in view \"%s\"): %lu\n\n",
    argv[2],
    (unsigned long)
    Libnucnet__Nuc__getNumberOfSpecies(
      Libnucnet__NucView__getNuc( p_view )
    )
  );

  /*============================================================================
  // Clean up.  Done.
  //==========================================================================*/

  Libnucnet__Nuc__free( p_my_nuclei );

  Libnucnet__NucView__free( p_view );

  /*============================================================================
  // Done!
  //==========================================================================*/

  return EXIT_SUCCESS;

}

/*##############################################################################
// modify_source().
//############################################################################*/

int
modify_species( Libnucnet__Species * p_species, double * p_addend )
{

  Libnucnet__Species__updateMassExcess(
    p_species,
    Libnucnet__Species__getMassExcess( p_species ) + *p_addend
  );

  Libnucnet__Species__updateSource(
    p_species,
    "Modified mass excess"
  );

  return 1;

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
    "\n  Index  Z    A    Name    Mass Excess (MeV)    Spin   Data Source\n"
  );
  fprintf(
    stdout,
    "  _____ ___  ___  ______  ___________________   ____   ___________\n\n"
  );

  /*============================================================================
  // Print out species data.
  //==========================================================================*/

  Libnucnet__Nuc__iterateSpecies(
    self,
    (Libnucnet__Species__iterateFunction) print_species,
    NULL
  );

  /*============================================================================
  // Print out the number of species.
  //==========================================================================*/

  fprintf( stdout, "\nThe collection has a total of %lu species.\n\n",
    (unsigned long) Libnucnet__Nuc__getNumberOfSpecies( self )
  );

  return;

}

/*##############################################################################
// print_species()
//############################################################################*/

int
print_species(
  Libnucnet__Species *p_species,
  void *p_data
)
{

  if( p_data ) {
    fprintf( stderr, "No extra data for this routine.\n" );
    exit( EXIT_FAILURE );
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

