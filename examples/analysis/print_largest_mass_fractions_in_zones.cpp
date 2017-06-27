////////////////////////////////////////////////////////////////////////////////
// This file was originally written by Tianhong Yu.
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
////////////////////////////////////////////////////////////////////////////////
 
////////////////////////////////////////////////////////////////////////////////
//! \file
//! \brief Example code to print out the largest mass fractions in zones.
//!    Output is zone number and the species and their mass fractions.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <string>
#include <iostream>
#include <Libnucnet.h>
#include <Libnucnet__Nuc.h>
#include <vector>

#include "nnt/auxiliary.h"
#include "nnt/iter.h"

typedef struct {
  const char *sName;
  double dMassFraction;
} species_struct;

typedef struct {
  Libnucnet__Zone *pZone;
  species_struct **pSpeciesStruct;
} Work; 

int
create_species_vector(
  Libnucnet__Species *,
  void *
);

int
compare_mass_fraction(
  const void *, const void *
);

//##############################################################################
// main().
//##############################################################################

int
main( int argc, char **argv )

{

  Libnucnet *p_my_nucnet;
  size_t i, i_number;
  
  species_struct **p_species_struct;
  Work *p_work; 

  if( argc != 4 )
  {
    fprintf(
      stderr,
      "\nUsage: %s xml_file top_number zone_xpath\n\n",
      argv[0]
    );
    fprintf(
      stderr,
      "  xml_file = network xml output file\n\n"
    );
    fprintf(
      stderr,
      "  top_number = how many mass fractions to print out\n\n"
    );
    fprintf(
      stderr,
      "  zone_xpath = XPath expression to select zones\n\n"
    );
    return EXIT_FAILURE;
  }

  p_my_nucnet =
    Libnucnet__new_from_xml(
      argv[1],
      NULL,
      NULL,
      argv[3]
    );

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function)
      nnt::zone_compare_by_first_label
  );

  i_number = 
    Libnucnet__Nuc__getNumberOfSpecies(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) )
    );

  p_species_struct = 
    ( species_struct ** )
    malloc(
      sizeof( species_struct * ) * i_number
    );

  for( i = 0; i < i_number; i++ ) {

    p_species_struct[i] = 
      ( species_struct * ) malloc( sizeof( species_struct ) );

  }

  p_work = ( Work * ) malloc( sizeof( Work ) );

  if( !p_species_struct || !p_work )
    LIBNUCNET__ERROR( "Couldn't allocate memory!" );

  p_work->pSpeciesStruct = p_species_struct;

  //============================================================================
  // Iterate the zones.
  //============================================================================

  BOOST_FOREACH( nnt::Zone zone, nnt::make_zone_list( p_my_nucnet ) )
  {

    std::cout << Libnucnet__Zone__getLabel( zone.getNucnetZone(), 1 );

    p_work->pZone = zone.getNucnetZone();

    Libnucnet__Nuc__iterateSpecies(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
      (Libnucnet__Species__iterateFunction) create_species_vector,
      p_work
    ); 

    qsort( 
      p_species_struct, 
      Libnucnet__Nuc__getNumberOfSpecies(
        Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) )
      ),
      sizeof( species_struct * ),
      compare_mass_fraction
    );

    std::cout << " ";

    for( i = 0; i < (size_t) atoi( argv[2] ); i++ ) {
      std::cout << p_work->pSpeciesStruct[i]->sName << " ";
      std::cout << p_work->pSpeciesStruct[i]->dMassFraction << " ";
    }

    std::cout << std::endl;

  }

  //============================================================================
  // Clean up.
  //============================================================================

  for( i = 0; i < i_number; i++ ) 
    free( p_species_struct[i] ); 

  free( p_species_struct );
  free( p_work );

  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}

int
create_species_vector(
  Libnucnet__Species *p_species,
  void *p_data
)
{

  Work *p_work = (Work *) p_data;

  size_t i = Libnucnet__Species__getIndex( p_species );

  p_work->pSpeciesStruct[i]->sName = Libnucnet__Species__getName( p_species );
  p_work->pSpeciesStruct[i]->dMassFraction = 
    Libnucnet__Zone__getSpeciesAbundance(
      p_work->pZone,
      p_species
    ) *
    (double) Libnucnet__Species__getA( p_species ); 

  return 1;

}

int
compare_mass_fraction(
  const void *a, const void *b
)
{

  species_struct *p_a = *(species_struct * const * ) a;
  species_struct *p_b = *(species_struct * const * ) b;

  if(
    p_a->dMassFraction > p_b->dMassFraction 
  ) {
    return -1;
  } else if( p_a->dMassFraction < p_b->dMassFraction ) {
    return 1;
  } else {
    return 0;
  }

}
    
