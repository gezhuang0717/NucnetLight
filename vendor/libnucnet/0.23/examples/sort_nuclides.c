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
//       xml file, sort the species according to a user-supplied data compare
//       function, print out the data about the species,
//       and clear the structure and free the allocated memory.
//     </abstract>
//   </description>
// </file>
//////////////////////////////////////////////////////////////////////////////*/
 
#include <Libnucnet__Nuc.h>

void print_nuclei( Libnucnet__Nuc * );

int print_species( Libnucnet__Species *, void * );

int
reverse_sort_function(
  const Libnucnet__Species *,
  const Libnucnet__Species *
);

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
  // Do a reverse sort on the nuclides.
  //==========================================================================*/

  Libnucnet__Nuc__setSpeciesCompareFunction(
    p_my_nuclei, 
    (Libnucnet__Species__compare_function) reverse_sort_function
  );

  Libnucnet__Nuc__sortSpecies( p_my_nuclei );

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

/*##############################################################################
// reverse_sort_function()
//############################################################################*/

int
reverse_sort_function(
  const Libnucnet__Species *p_species1,
  const Libnucnet__Species *p_species2
)
{

  int i;

  if( 
      Libnucnet__Species__getZ( p_species1 ) >
      Libnucnet__Species__getZ( p_species2 )
  )
    return -1;
  else if(
      Libnucnet__Species__getZ( p_species1 ) <
      Libnucnet__Species__getZ( p_species2 )
  )
    return 1;

  if( 
      Libnucnet__Species__getA( p_species1 ) >
      Libnucnet__Species__getA( p_species2 )
  )
    return -1;
  else if(
      Libnucnet__Species__getA( p_species1 ) <
      Libnucnet__Species__getA( p_species2 )
  )
    return 1;

  i =
    strcmp(
      Libnucnet__Species__getName( p_species1 ),
      Libnucnet__Species__getName( p_species2 )
    );

  if( i == 0 ) {
    return 0;
  } else {
    return -GSL_SIGN( i );
  }

}
   
